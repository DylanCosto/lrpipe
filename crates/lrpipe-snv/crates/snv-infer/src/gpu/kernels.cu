// Kernels of the Rust GPU backend, compiled at run time by NVRTC (no CUDA headers needed).
//
// Matrices are row-major `float`. Every kernel runs 256 threads per block.

#define SELU_ALPHA 1.6732632423543772f
#define SELU_SCALE 1.0507009873554805f

// Activations after a product's bias (and residual).
#define ACT_NONE 0
#define ACT_RELU 1
#define ACT_SELU 2

#define NEG_INF __int_as_float(0xff800000)

// x rounded to TF32 (10 mantissa bits, to nearest, ties away), as tensor cores take float32
// operands; products of two such values are exact in float32.
__device__ __forceinline__ float tf32(float x) {
#if __CUDA_ARCH__ >= 800
    unsigned r;
    asm("cvt.rna.tf32.f32 %0, %1;" : "=r"(r) : "f"(x));
    return __uint_as_float(r);
#else
    // Before Ampere, which has no such instruction, the same rounding in integer arithmetic: half a
    // unit of the 13 dropped mantissa bits added to the magnitude, then the bits dropped (a carry
    // reaches the exponent as it should). Infinities and NaNs stay as they are.
    unsigned u = __float_as_uint(x);
    if ((u & 0x7f800000u) != 0x7f800000u) u = (u + 0x1000u) & 0xffffe000u;
    return __uint_as_float(u);
#endif
}

__device__ __forceinline__ float activate(float x, int act) {
    if (act == ACT_RELU) return fmaxf(x, 0.0f);
    if (act == ACT_SELU) return x > 0.0f ? SELU_SCALE * x : SELU_SCALE * SELU_ALPHA * expm1f(x);
    return x;
}

// Where the left operand of a product comes from.
struct Lhs {
    const float* a;
    // Dense: element (row, k) is a[row * lda + k].
    int lda;
    // Convolution (implicit im2col over a [image][h][w][c] input): row is (image, oy, ox) and
    // k is (ky, kx, c) of a 3x3 kernel with padding 1.
    int h, w, c, oh, ow, stride;
};

// One thread's share of loading a BM x 8 tile of the left operand: rows tid / 8 + 32 i at
// k = k0 + tid % 8. What depends on the row alone is worked out once, before the k loop.
template <int BM, bool CONV>
struct LhsLoader {
    static constexpr int ROWS = BM * 8 / 256;
    const float* base[ROWS];
    // Convolution: input position of the kernel's top-left tap (may be -1).
    int y0[ROWS], x0[ROWS];
    bool valid[ROWS];
    int kk;
    // Convolution: tap (ky * 3 + kx) and channel of k0 + kk.
    int tap, ch;

    __device__ __forceinline__ LhsLoader(const Lhs& l, int row0, int m) {
        kk = threadIdx.x % 8;
#pragma unroll
        for (int i = 0; i < ROWS; i++) {
            int row = row0 + threadIdx.x / 8 + 32 * i;
            valid[i] = row < m;
            if (!valid[i]) row = 0;
            if (CONV) {
                int image = row / (l.oh * l.ow), pos = row - image * (l.oh * l.ow);
                int oy = pos / l.ow, ox = pos - oy * l.ow;
                base[i] = l.a + (size_t)image * l.h * l.w * l.c;
                y0[i] = oy * l.stride - 1;
                x0[i] = ox * l.stride - 1;
            } else {
                base[i] = l.a + (size_t)row * l.lda;
                y0[i] = x0[i] = 0;
            }
        }
        tap = CONV ? kk / l.c : 0;
        ch = CONV ? kk - tap * l.c : 0;
    }

    // The values at k = k0 + kk of this thread's rows, then moves on by 8.
    __device__ __forceinline__ void fetch(float (&v)[ROWS], const Lhs& l, int k0, int kdim) {
        const int k = k0 + kk;
        const int ky = tap / 3, kx = tap - 3 * (tap / 3);
#pragma unroll
        for (int i = 0; i < ROWS; i++) {
            v[i] = 0.0f;
            if (valid[i] && k < kdim) {
                if (CONV) {
                    int iy = y0[i] + ky, ix = x0[i] + kx;
                    if (iy >= 0 && ix >= 0 && iy < l.h && ix < l.w) v[i] = base[i][((size_t)iy * l.w + ix) * l.c + ch];
                } else {
                    v[i] = base[i][k];
                }
            }
        }
        if (CONV) {
            ch += 8;
            while (ch >= l.c) {
                ch -= l.c;
                tap++;
            }
        }
    }
};

__device__ __forceinline__ float epilogue(float acc, int row, int col, int ldc, const float* scale, const float* bias,
                                          const float* residual, int act);

// Split-K: block z of the grid sums k in [*begin, *end), in whole tiles of `tile`.
__device__ __forceinline__ void split_range(int k, int tile, int* begin, int* end) {
    int chunk = (k + gridDim.z - 1) / gridDim.z;
    chunk = (chunk + tile - 1) / tile * tile;
    *begin = blockIdx.z * chunk;
    *end = min(k, *begin + chunk);
}

// Stores a sum: with split-K (gridDim.z > 1) raw into slice z of `out` (m x ldc each), for
// splitk_reduce; otherwise through the epilogue.
__device__ __forceinline__ void store_sum(float* out, float acc, int row, int col, int m, int ldc, const float* scale,
                                          const float* bias, const float* residual, int act) {
    if (gridDim.z > 1) out[((size_t)blockIdx.z * m + row) * ldc + col] = acc;
    else out[(size_t)row * ldc + col] = epilogue(acc, row, col, ldc, scale, bias, residual, act);
}

// act(acc * scale + bias (+ residual)); `scale`, a batch norm's, may be null.
__device__ __forceinline__ float epilogue(float acc, int row, int col, int ldc, const float* scale, const float* bias,
                                          const float* residual, int act) {
    float v = (scale ? acc * scale[col] : acc) + bias[col];
    if (residual) v += residual[(size_t)row * ldc + col];
    return activate(v, act);
}

// out[m x n] = act(lhs[m x k] * b[k x n] + bias (+ residual)), b dense with row length n.
// Block tile BM x BN, 16 x 16 threads, each owning rows ty + 16 i and columns tx + 16 j.
template <int BM, int BN, bool CONV>
__device__ __forceinline__ void gemm(const Lhs& lhs, const float* __restrict__ b, float* __restrict__ out,
                                     int ldc, const float* __restrict__ scale, const float* __restrict__ bias, const float* __restrict__ residual,
                                     int m, int n, int k, int act, int round) {
    constexpr int BK = 8, TM = BM / 16, TN = BN / 16;
    __shared__ float as[BK][BM + 4];
    __shared__ float bs[BK][BN];
    const int tid = threadIdx.x, tx = tid % 16, ty = tid / 16;
    const int row0 = blockIdx.y * BM, col0 = blockIdx.x * BN;
    float acc[TM][TN];
#pragma unroll
    for (int i = 0; i < TM; i++)
#pragma unroll
        for (int j = 0; j < TN; j++) acc[i][j] = 0.0f;
    LhsLoader<BM, CONV> loader(lhs, row0, m);
    int kbegin, kend;
    split_range(k, BK, &kbegin, &kend);
    for (int k0 = kbegin; k0 < kend; k0 += BK) {
        {
            float v[LhsLoader<BM, CONV>::ROWS];
            loader.fetch(v, lhs, k0, kend);
#pragma unroll
            for (int i = 0; i < LhsLoader<BM, CONV>::ROWS; i++) as[tid % 8][tid / 8 + 32 * i] = round ? tf32(v[i]) : v[i];
        }
#pragma unroll
        for (int e = tid; e < BK * BN; e += 256) {
            int kk = e / BN, col = e % BN;
            float v = (k0 + kk < kend && col0 + col < n) ? b[(size_t)(k0 + kk) * n + col0 + col] : 0.0f;
            bs[kk][col] = round ? tf32(v) : v;
        }
        __syncthreads();
#pragma unroll
        for (int kk = 0; kk < BK; kk++) {
            float av[TM], bv[TN];
#pragma unroll
            for (int i = 0; i < TM; i++) av[i] = as[kk][ty + 16 * i];
#pragma unroll
            for (int j = 0; j < TN; j++) bv[j] = bs[kk][tx + 16 * j];
#pragma unroll
            for (int i = 0; i < TM; i++)
#pragma unroll
                for (int j = 0; j < TN; j++) acc[i][j] = fmaf(av[i], bv[j], acc[i][j]);
        }
        __syncthreads();
    }
#pragma unroll
    for (int i = 0; i < TM; i++) {
        int row = row0 + ty + 16 * i;
        if (row >= m) continue;
#pragma unroll
        for (int j = 0; j < TN; j++) {
            int col = col0 + tx + 16 * j;
            if (col >= n) continue;
            store_sum(out, acc[i][j], row, col, m, ldc, scale, bias, residual, act);
        }
    }
}

// As gemm for BM = 128, faster: each thread owns rows 4 ty..4 ty+3 and 64 + the same, and
// columns 4 tx..4 tx+3 (and 64 + those when BN = 128), read from shared memory as float4; the
// next tile is fetched into registers while this one is multiplied (two shared buffers).
template <int BN, bool CONV>
__device__ __forceinline__ void gemm128(const Lhs& lhs, const float* __restrict__ b, float* __restrict__ out,
                                        int ldc, const float* __restrict__ scale, const float* __restrict__ bias, const float* __restrict__ residual,
                                        int m, int n, int k, int act, int round) {
    constexpr int BM = 128, BK = 8, CN = BN / 64, BLOADS = BK * BN / 256;
    __shared__ __align__(16) float as[2][BK][BM + 4];
    __shared__ __align__(16) float bs[2][BK][BN];
    const int tid = threadIdx.x, tx = tid % 16, ty = tid / 16;
    const int row0 = blockIdx.y * BM, col0 = blockIdx.x * BN;
    float acc[8][4 * CN];
#pragma unroll
    for (int i = 0; i < 8; i++)
#pragma unroll
        for (int j = 0; j < 4 * CN; j++) acc[i][j] = 0.0f;
    LhsLoader<BM, CONV> loader(lhs, row0, m);
    float av[LhsLoader<BM, CONV>::ROWS], bv[BLOADS];
    int kbegin, kend;
    split_range(k, BK, &kbegin, &kend);
    auto fetch = [&](int k0) {
        loader.fetch(av, lhs, k0, kend);
#pragma unroll
        for (int j = 0; j < BLOADS; j++) {
            int e = tid + 256 * j, kk = e / BN, col = e % BN;
            bv[j] = (k0 + kk < kend && col0 + col < n) ? b[(size_t)(k0 + kk) * n + col0 + col] : 0.0f;
        }
    };
    auto store = [&](int buffer) {
#pragma unroll
        for (int i = 0; i < LhsLoader<BM, CONV>::ROWS; i++) as[buffer][tid % 8][tid / 8 + 32 * i] = round ? tf32(av[i]) : av[i];
#pragma unroll
        for (int j = 0; j < BLOADS; j++) {
            int e = tid + 256 * j;
            bs[buffer][e / BN][e % BN] = round ? tf32(bv[j]) : bv[j];
        }
    };
    fetch(kbegin);
    store(0);
    __syncthreads();
    int buffer = 0;
    for (int k0 = kbegin; k0 < kend; k0 += BK) {
        const bool more = k0 + BK < kend;
        if (more) fetch(k0 + BK);
#pragma unroll
        for (int kk = 0; kk < BK; kk++) {
            float a[8], c[4 * CN];
            *(float4*)&a[0] = *(const float4*)&as[buffer][kk][4 * ty];
            *(float4*)&a[4] = *(const float4*)&as[buffer][kk][64 + 4 * ty];
#pragma unroll
            for (int j = 0; j < CN; j++) *(float4*)&c[4 * j] = *(const float4*)&bs[buffer][kk][64 * j + 4 * tx];
#pragma unroll
            for (int i = 0; i < 8; i++)
#pragma unroll
                for (int j = 0; j < 4 * CN; j++) acc[i][j] = fmaf(a[i], c[j], acc[i][j]);
        }
        if (more) {
            store(buffer ^ 1);
            __syncthreads();
            buffer ^= 1;
        }
    }
#pragma unroll
    for (int i = 0; i < 8; i++) {
        int row = row0 + (i < 4 ? 4 * ty + i : 64 + 4 * ty + i - 4);
        if (row >= m) continue;
#pragma unroll
        for (int j = 0; j < 4 * CN; j++) {
            int col = col0 + 64 * (j / 4) + 4 * tx + j % 4;
            if (col < n) store_sum(out, acc[i][j], row, col, m, ldc, scale, bias, residual, act);
        }
    }
}

#define DENSE_GEMM(NAME, CALL)                                                                           \
    extern "C" __global__ void __launch_bounds__(256)                                                     \
        NAME(const float* a, int lda, const float* b, float* out, int ldc, const float* bias,             \
             const float* residual, int m, int n, int k, int act, int round) {                            \
        Lhs lhs = {a, lda, 0, 0, 0, 0, 0, 0};                                                             \
        CALL(lhs, b, out, ldc, nullptr, bias, residual, m, n, k, act, round);                                      \
    }

#define CONV_GEMM(NAME, BM, BN)                                                                          \
    extern "C" __global__ void __launch_bounds__(256)                                                     \
        NAME(const float* input, int h, int w, int c, int oh, int ow, int stride, const float* b,         \
             float* out, const float* scale, const float* bias, const float* residual, int m, int n,      \
             int act, int round) {                                                                        \
        Lhs lhs = {input, 0, h, w, c, oh, ow, stride};                                                    \
        gemm128<BN, true>(lhs, b, out, n, scale, bias, residual, m, n, 9 * c, act, round);                \
    }

DENSE_GEMM(gemm_128x128, (gemm128<128, false>))
DENSE_GEMM(gemm_128x64, (gemm128<64, false>))
DENSE_GEMM(gemm_64x32, (gemm<64, 32, false>))
CONV_GEMM(conv_128x128, 128, 128)
CONV_GEMM(conv_128x64, 128, 64)

// ---- TF32 tensor cores (Ampere, 8.0, and later): the products libtorch's cuDNN runs on them. ----

// d += a * b for one m16n8k8 tile of TF32 operands (float bits, already rounded by tf32()), with
// float32 sums. Before Ampere the host never launches the kernels that use it.
__device__ __forceinline__ void mma_tf32(float (&d)[4], const unsigned (&a)[4], const unsigned (&b)[2]) {
#if __CUDA_ARCH__ >= 800
    asm volatile(
        "mma.sync.aligned.m16n8k8.row.col.f32.tf32.tf32.f32 {%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, "
        "{%0,%1,%2,%3};"
        : "+f"(d[0]), "+f"(d[1]), "+f"(d[2]), "+f"(d[3])
        : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]), "r"(b[0]), "r"(b[1]));
#else
    __trap();
#endif
}

// As gemm128 with TF32 operands, on the tensor cores: block tile 128 x BN x 16, eight warps as
// 2 x 4, each owning 64 rows and BN / 4 columns. The next tile is fetched into registers while
// this one is multiplied (two shared buffers), and the paddings keep the fragment reads free of
// shared-memory bank conflicts.
template <int BN, bool CONV>
__device__ __forceinline__ void gemm_tf32(const Lhs& lhs, const float* __restrict__ b, float* __restrict__ out,
                                         int ldc, const float* __restrict__ scale, const float* __restrict__ bias,
                                         const float* __restrict__ residual, int m, int n, int k, int act) {
    constexpr int BM = 128, BK = 16, WN = BN / 4, MT = 4, NT = WN / 8;
    constexpr int AROWS = LhsLoader<BM, CONV>::ROWS, BLOADS = BK * BN / 256;
    __shared__ __align__(16) float as[2][BM][BK + 4];
    __shared__ __align__(16) float bs[2][BK][BN + 8];
    const int tid = threadIdx.x, lane = tid % 32, warp = tid / 32;
    const int wm = warp / 4, wn = warp % 4, group = lane / 4, tig = lane % 4;
    const int row0 = blockIdx.y * BM, col0 = blockIdx.x * BN;
    float acc[MT][NT][4];
#pragma unroll
    for (int i = 0; i < MT; i++)
#pragma unroll
        for (int j = 0; j < NT; j++)
#pragma unroll
            for (int q = 0; q < 4; q++) acc[i][j][q] = 0.0f;
    LhsLoader<BM, CONV> loader(lhs, row0, m);
    float av[BK / 8][AROWS], bv[BLOADS];
    int kbegin, kend;
    split_range(k, BK, &kbegin, &kend);
    auto fetch = [&](int k0) {
#pragma unroll
        for (int j = 0; j < BK / 8; j++) loader.fetch(av[j], lhs, k0 + 8 * j, kend);
#pragma unroll
        for (int j = 0; j < BLOADS; j++) {
            int e = tid + 256 * j, kk = e / BN, col = col0 + e % BN;
            bv[j] = (k0 + kk < kend && col < n) ? b[(size_t)(k0 + kk) * n + col] : 0.0f;
        }
    };
    auto store = [&](int buffer) {
#pragma unroll
        for (int j = 0; j < BK / 8; j++)
#pragma unroll
            for (int i = 0; i < AROWS; i++) as[buffer][tid / 8 + 32 * i][8 * j + tid % 8] = tf32(av[j][i]);
#pragma unroll
        for (int j = 0; j < BLOADS; j++) {
            int e = tid + 256 * j;
            bs[buffer][e / BN][e % BN] = tf32(bv[j]);
        }
    };
    fetch(kbegin);
    store(0);
    __syncthreads();
    int buffer = 0;
    for (int k0 = kbegin; k0 < kend; k0 += BK) {
        const bool more = k0 + BK < kend;
        if (more) fetch(k0 + BK);
#pragma unroll
        for (int ks = 0; ks < BK; ks += 8) {
            unsigned a[MT][4], bf[NT][2];
#pragma unroll
            for (int i = 0; i < MT; i++) {
                const int r = wm * 64 + i * 16 + group;
                a[i][0] = __float_as_uint(as[buffer][r][ks + tig]);
                a[i][1] = __float_as_uint(as[buffer][r + 8][ks + tig]);
                a[i][2] = __float_as_uint(as[buffer][r][ks + tig + 4]);
                a[i][3] = __float_as_uint(as[buffer][r + 8][ks + tig + 4]);
            }
#pragma unroll
            for (int j = 0; j < NT; j++) {
                const int cn = wn * WN + j * 8 + group;
                bf[j][0] = __float_as_uint(bs[buffer][ks + tig][cn]);
                bf[j][1] = __float_as_uint(bs[buffer][ks + tig + 4][cn]);
            }
#pragma unroll
            for (int i = 0; i < MT; i++)
#pragma unroll
                for (int j = 0; j < NT; j++) mma_tf32(acc[i][j], a[i], bf[j]);
        }
        if (more) {
            store(buffer ^ 1);
            __syncthreads();
            buffer ^= 1;
        }
    }
#pragma unroll
    for (int i = 0; i < MT; i++)
#pragma unroll
        for (int j = 0; j < NT; j++)
#pragma unroll
            for (int q = 0; q < 4; q++) {
                int row = row0 + wm * 64 + i * 16 + group + (q / 2) * 8;
                int col = col0 + wn * WN + j * 8 + tig * 2 + q % 2;
                if (row < m && col < n) store_sum(out, acc[i][j][q], row, col, m, ldc, scale, bias, residual, act);
            }
}

// The same parameters as DENSE_GEMM and CONV_GEMM (`round` is implied).
#define DENSE_TF32(NAME, BN)                                                                             \
    extern "C" __global__ void __launch_bounds__(256)                                                     \
        NAME(const float* a, int lda, const float* b, float* out, int ldc, const float* bias,             \
             const float* residual, int m, int n, int k, int act, int round) {                            \
        Lhs lhs = {a, lda, 0, 0, 0, 0, 0, 0};                                                             \
        gemm_tf32<BN, false>(lhs, b, out, ldc, nullptr, bias, residual, m, n, k, act);                    \
    }

#define CONV_TF32(NAME, BN)                                                                              \
    extern "C" __global__ void __launch_bounds__(256)                                                     \
        NAME(const float* input, int h, int w, int c, int oh, int ow, int stride, const float* b,         \
             float* out, const float* scale, const float* bias, const float* residual, int m, int n,      \
             int act, int round) {                                                                        \
        Lhs lhs = {input, 0, h, w, c, oh, ow, stride};                                                    \
        gemm_tf32<BN, true>(lhs, b, out, n, scale, bias, residual, m, n, 9 * c, act);                     \
    }

DENSE_TF32(gemm_tf32_128x128, 128)
DENSE_TF32(gemm_tf32_128x64, 64)
CONV_TF32(conv_tf32_128x128, 128)
CONV_TF32(conv_tf32_128x64, 64)

__device__ __forceinline__ float sigmoid(float x) { return 1.0f / (1.0f + expf(-x)); }

// One step of both directions (blockIdx.z) of a bidirectional LSTM layer with hidden size h over
// `steps` steps: step s of the forward direction is position s, of the reverse one steps - 1 - s.
//   pre:  [n][steps][2][4][h] input projections plus both biases (PyTorch's gate order i, f, g, o)
//   whh:  [2][h][4 h], W_hh transposed per direction
//   out:  [n][steps][2][h] hidden states; position `prev` is read, `t` written
//   cell: [2][n][h]
// A block owns 64 rows and 32 hidden units (all four gates of each), a thread 4 adjacent rows
// and 2 adjacent units, so each step of the sum reads its operands from shared memory as one
// float4 and four float2. Each sum runs over k in order, as before the tiling changed.
extern "C" __global__ void __launch_bounds__(256)
    lstm_step(const float* __restrict__ pre, const float* __restrict__ whh, float* out, float* __restrict__ cell,
              int n, int steps, int h, int s, int round) {
    constexpr int RB = 64, UB = 32, KT = 16;
    __shared__ __align__(16) float hs[KT][RB];
    __shared__ __align__(16) float ws[KT][4][UB];
    const int d = blockIdx.z;
    const int t = d == 0 ? s : steps - 1 - s;
    const int prev = d == 0 ? t - 1 : t + 1;
    const int tid = threadIdx.x, tx = tid % 16, ty = tid / 16;
    const int r0 = blockIdx.x * RB, u0 = blockIdx.y * UB;
    float acc[4][2][4];
#pragma unroll
    for (int i = 0; i < 4; i++)
#pragma unroll
        for (int j = 0; j < 2; j++)
#pragma unroll
            for (int g = 0; g < 4; g++) acc[i][j][g] = 0.0f;
    if (s > 0) {
        const float* w = whh + (size_t)d * h * 4 * h;
        for (int k0 = 0; k0 < h; k0 += KT) {
            for (int e = tid; e < RB * KT; e += 256) {
                int r = e / KT, kk = e % KT, row = r0 + r;
                float v = (row < n && k0 + kk < h)
                              ? out[((size_t)row * steps + prev) * 2 * h + (size_t)d * h + k0 + kk]
                              : 0.0f;
                hs[kk][r] = round ? tf32(v) : v;
            }
            for (int e = tid; e < KT * 4 * UB; e += 256) {
                int kk = e / (4 * UB), g = (e / UB) % 4, u = e % UB;
                float v = (k0 + kk < h && u0 + u < h) ? w[(size_t)(k0 + kk) * 4 * h + g * h + u0 + u] : 0.0f;
                ws[kk][g][u] = round ? tf32(v) : v;
            }
            __syncthreads();
#pragma unroll
            for (int kk = 0; kk < KT; kk++) {
                const float4 hv = *(const float4*)&hs[kk][4 * ty];
                float2 wv[4];
#pragma unroll
                for (int g = 0; g < 4; g++) wv[g] = *(const float2*)&ws[kk][g][2 * tx];
                const float hr[4] = {hv.x, hv.y, hv.z, hv.w};
#pragma unroll
                for (int i = 0; i < 4; i++)
#pragma unroll
                    for (int g = 0; g < 4; g++) {
                        acc[i][0][g] = fmaf(hr[i], wv[g].x, acc[i][0][g]);
                        acc[i][1][g] = fmaf(hr[i], wv[g].y, acc[i][1][g]);
                    }
            }
            __syncthreads();
        }
    }
#pragma unroll
    for (int i = 0; i < 4; i++) {
        int row = r0 + 4 * ty + i;
        if (row >= n) continue;
#pragma unroll
        for (int j = 0; j < 2; j++) {
            int u = u0 + 2 * tx + j;
            if (u >= h) continue;
            const float* p = pre + (((size_t)row * steps + t) * 2 + d) * 4 * h;
            float* c = cell + ((size_t)d * n + row) * h + u;
            float gi = sigmoid(acc[i][j][0] + p[u]);
            float gf = sigmoid(acc[i][j][1] + p[h + u]);
            float gg = tanhf(acc[i][j][2] + p[2 * h + u]);
            float go = sigmoid(acc[i][j][3] + p[3 * h + u]);
            float cv = (s > 0 ? gf * *c : 0.0f) + gi * gg;
            *c = cv;
            out[((size_t)row * steps + t) * 2 * h + (size_t)d * h + u] = go * tanhf(cv);
        }
    }
}

// x[i] = (float)in[i] / divisor.
extern "C" __global__ void int8_to_float(const signed char* in, float* x, long long len, float divisor) {
    for (long long i = blockIdx.x * 256LL + threadIdx.x; i < len; i += (long long)gridDim.x * 256) x[i] = in[i] / divisor;
}

extern "C" __global__ void int32_to_float(const int* in, float* x, long long len) {
    for (long long i = blockIdx.x * 256LL + threadIdx.x; i < len; i += (long long)gridDim.x * 256) x[i] = (float)in[i];
}

// PyramidPolling: TF-style "same" max pools of each [h][w][c] map to 3x3, 2x2 and 1x1 grids, the
// extra padding row or column at the end and counting as zero, each flattened [row][col][c].
// out: [images][14 c].
extern "C" __global__ void pyramid_pool(const float* in, float* out, int images, int h, int w, int c) {
    const int features = 14 * c;
    for (long long e = blockIdx.x * 256LL + threadIdx.x; e < (long long)images * features;
         e += (long long)gridDim.x * 256) {
        int image = e / features, f = e % features;
        int size, cell;
        if (f < 9 * c) { size = 3; cell = f / c; }
        else if (f < 13 * c) { size = 2; cell = (f - 9 * c) / c; }
        else { size = 1; cell = 0; }
        int ch = f % c;
        int wh = (h + size - 1) / size, ww = (w + size - 1) / size;
        int out_h = (h + wh - 1) / wh, out_w = (w + ww - 1) / ww;
        int pad_h = out_h * wh - h, pad_w = out_w * ww - w;
        pad_h = pad_h > 0 ? pad_h : 0;
        pad_w = pad_w > 0 ? pad_w : 0;
        int oy = cell / out_w, ox = cell % out_w;
        int y0 = oy * wh - pad_h / 2, x0 = ox * ww - pad_w / 2;
        int y1 = y0 + wh, x1 = x0 + ww;
        bool padded = y0 < 0 || x0 < 0 || y1 > h || x1 > w;
        y0 = y0 < 0 ? 0 : y0; x0 = x0 < 0 ? 0 : x0;
        y1 = y1 > h ? h : y1; x1 = x1 > w ? w : x1;
        float m = padded ? 0.0f : NEG_INF;
        const float* map = in + (size_t)image * h * w * c;
        for (int y = y0; y < y1; y++)
            for (int x = x0; x < x1; x++) m = fmaxf(m, map[((size_t)y * w + x) * c + ch]);
        out[(size_t)image * features + f] = m;
    }
}

// Softmax over each row of logits [n][classes], written to out[row * ldo + offset + class].
extern "C" __global__ void softmax_rows(const float* logits, float* out, int n, int classes, int ldo, int offset) {
    int row = blockIdx.x * 256 + threadIdx.x;
    if (row >= n) return;
    const float* l = logits + (size_t)row * classes;
    float m = NEG_INF;
    for (int i = 0; i < classes; i++) m = fmaxf(m, l[i]);
    float sum = 0.0f;
    for (int i = 0; i < classes; i++) sum += expf(l[i] - m);
    for (int i = 0; i < classes; i++) out[(size_t)row * ldo + offset + i] = expf(l[i] - m) / sum;
}

// out = act(sum of the `splits` partial sums in ws (each m x n) + bias), for split-K products.
extern "C" __global__ void splitk_reduce(const float* ws, float* out, int m, int n, int ldc, int splits,
                                         const float* bias, int act) {
    for (long long e = blockIdx.x * 256LL + threadIdx.x; e < (long long)m * n; e += (long long)gridDim.x * 256) {
        int row = e / n, col = e % n;
        float v = 0.0f;
        for (int z = 0; z < splits; z++) v += ws[((size_t)z * m + row) * n + col];
        out[(size_t)row * ldc + col] = activate(v + bias[col], act);
    }
}
