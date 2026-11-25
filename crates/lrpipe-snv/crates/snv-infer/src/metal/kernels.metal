// Kernels of the Rust Metal backend, compiled at run time by the Metal framework.
//
// Matrices are row-major `float`, in float32 throughout (Apple GPUs have no TF32), with precise
// math: the library is compiled with fast math off.

#include <metal_stdlib>
using namespace metal;

// Activations after a product's bias (and residual).
#define ACT_NONE 0
#define ACT_RELU 1
#define ACT_SELU 2

constant float SELU_ALPHA = 1.6732632423543772f;
constant float SELU_SCALE = 1.0507009873554805f;

// exp(x) - 1, accurate near 0 (Kahan's: the error of exp(x) cancels in (u - 1) / log(u)).
// Metal's library has no expm1.
inline float expm1_(float x) {
    float u = precise::exp(x);
    if (u == 1.0f) return x;
    float um1 = u - 1.0f;
    if (um1 == -1.0f) return -1.0f;
    return um1 * x / precise::log(u);
}

inline float activate(float x, int act) {
    if (act == ACT_RELU) return max(x, 0.0f);
    if (act == ACT_SELU) return x > 0.0f ? SELU_SCALE * x : SELU_SCALE * SELU_ALPHA * expm1_(x);
    return x;
}

// The parameters of the products, as `metal::ProductParams` lays them out.
struct ProductParams {
    int m, n, k;
    // Rows of the dense left operand, the output and the residual are this many values apart.
    int lda, ldc, ldr;
    int act;
    int has_bias, has_residual;
    // Convolution (implicit im2col over a [image][h][w][c] input): row is (image, oy, ox) and k is
    // (ky, kx, c) of a 3x3 kernel with padding 1.
    int h, w, c, oh, ow, stride;
};

// One thread's share of loading a BM x 8 tile of the left operand: rows tid / 8 + (THREADS / 8) i
// at k = k0 + tid % 8. What depends on the row alone is worked out once, before the k loop.
template <int BM, int THREADS, bool CONV>
struct LhsLoader {
    static constant constexpr int ROWS = BM * 8 / THREADS;
    device const float* base[ROWS];
    // Convolution: input position of the kernel's top-left tap (may be -1).
    int y0[ROWS], x0[ROWS];
    bool valid[ROWS];
    int kk;
    // Convolution: tap (ky * 3 + kx) and channel of k0 + kk.
    int tap, ch;

    LhsLoader(device const float* a, constant ProductParams& p, int row0, int k0, uint tid) {
        kk = tid % 8;
        for (int i = 0; i < ROWS; i++) {
            int row = row0 + tid / 8 + (THREADS / 8) * i;
            valid[i] = row < p.m;
            if (!valid[i]) row = 0;
            if (CONV) {
                int image = row / (p.oh * p.ow), pos = row - image * (p.oh * p.ow);
                int oy = pos / p.ow, ox = pos - oy * p.ow;
                base[i] = a + (size_t)image * p.h * p.w * p.c;
                y0[i] = oy * p.stride - 1;
                x0[i] = ox * p.stride - 1;
            } else {
                base[i] = a + (size_t)row * p.lda;
                y0[i] = x0[i] = 0;
            }
        }
        tap = CONV ? (k0 + kk) / p.c : 0;
        ch = CONV ? (k0 + kk) - tap * p.c : 0;
    }

    // The values at k = k0 + kk of this thread's rows, then moves on by 8.
    void fetch(thread float* v, constant ProductParams& p, int k0, int kend) {
        const int k = k0 + kk;
        const int ky = tap / 3, kx = tap - 3 * (tap / 3);
        for (int i = 0; i < ROWS; i++) {
            v[i] = 0.0f;
            if (valid[i] && k < kend) {
                if (CONV) {
                    int iy = y0[i] + ky, ix = x0[i] + kx;
                    if (iy >= 0 && ix >= 0 && iy < p.h && ix < p.w) v[i] = base[i][((size_t)iy * p.w + ix) * p.c + ch];
                } else {
                    v[i] = base[i][k];
                }
            }
        }
        if (CONV) {
            ch += 8;
            while (ch >= p.c) {
                ch -= p.c;
                tap++;
            }
        }
    }
};

// Split-K: group z of the grid sums k in [*begin, *end), in whole tiles of `tile`.
inline void split_range(int k, int tile, uint3 group, uint3 groups, thread int& begin, thread int& end) {
    int chunk = (k + groups.z - 1) / groups.z;
    chunk = (chunk + tile - 1) / tile * tile;
    begin = group.z * chunk;
    end = min(k, begin + chunk);
}

// Stores a sum: with split-K (grid depth > 1) raw into slice z of `out` (m x n each), for
// splitk_reduce; otherwise act(acc + bias (+ residual)).
inline void store_sum(device float* out, float acc, int row, int col, constant ProductParams& p,
                      device const float* bias, device const float* residual, uint3 group, uint3 groups) {
    if (groups.z > 1) {
        out[((size_t)group.z * p.m + row) * p.n + col] = acc;
        return;
    }
    float v = acc;
    if (p.has_bias) v += bias[col];
    if (p.has_residual) v += residual[(size_t)row * p.ldr + col];
    out[(size_t)row * p.ldc + col] = activate(v, p.act);
}

// out[m x n] = act(lhs[m x k] * b[k x n] + bias (+ residual)), b dense with rows of n, for small
// products. Block tile 64 x 32, 16 x 16 threads, each owning rows ty + 16 i and columns tx + 16 j.
template <bool CONV>
inline void gemm_small(device const float* a, device const float* b, device float* out, device const float* bias,
                       device const float* residual, constant ProductParams& p, threadgroup float* as,
                       threadgroup float* bs, uint3 group, uint3 groups, uint tid) {
    constexpr int BM = 64, BN = 32, BK = 8, TM = BM / 16, TN = BN / 16;
    const int tx = tid % 16, ty = tid / 16;
    const int row0 = group.y * BM, col0 = group.x * BN;
    float acc[TM][TN];
    for (int i = 0; i < TM; i++)
        for (int j = 0; j < TN; j++) acc[i][j] = 0.0f;
    int kbegin, kend;
    split_range(p.k, BK, group, groups, kbegin, kend);
    LhsLoader<BM, 256, CONV> loader(a, p, row0, kbegin, tid);
    for (int k0 = kbegin; k0 < kend; k0 += BK) {
        float v[LhsLoader<BM, 256, CONV>::ROWS];
        loader.fetch(v, p, k0, kend);
        for (int i = 0; i < LhsLoader<BM, 256, CONV>::ROWS; i++) as[(tid % 8) * BM + tid / 8 + 32 * i] = v[i];
        for (int e = tid; e < BK * BN; e += 256) {
            int kk = e / BN, col = e % BN;
            bs[kk * BN + col] = (k0 + kk < kend && col0 + col < p.n) ? b[(size_t)(k0 + kk) * p.n + col0 + col] : 0.0f;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
        for (int kk = 0; kk < BK; kk++) {
            float av[TM], bv[TN];
            for (int i = 0; i < TM; i++) av[i] = as[kk * BM + ty + 16 * i];
            for (int j = 0; j < TN; j++) bv[j] = bs[kk * BN + tx + 16 * j];
            for (int i = 0; i < TM; i++)
                for (int j = 0; j < TN; j++) acc[i][j] = fma(av[i], bv[j], acc[i][j]);
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    for (int i = 0; i < TM; i++) {
        int row = row0 + ty + 16 * i;
        if (row >= p.m) continue;
        for (int j = 0; j < TN; j++) {
            int col = col0 + tx + 16 * j;
            if (col < p.n) store_sum(out, acc[i][j], row, col, p, bias, residual, group, groups);
        }
    }
}

// As gemm_small for BM = 128, faster: each thread owns rows 4 ty..4 ty+3 and 64 + the same, and
// columns 4 tx..4 tx+3 (and 64 + those when BN = 128), read from threadgroup memory as float4;
// the next tile is fetched into registers while this one is multiplied (two buffers). The host
// runs BN = 64: at 128, the 64 sums per thread ran a quarter as fast on an M4 Pro.
constant constexpr int BM = 128, BK = 8, AS = BM + 4;

template <int BN, bool CONV>
inline void gemm128(device const float* a, device const float* b, device float* out, device const float* bias,
                    device const float* residual, constant ProductParams& p, threadgroup float* as,
                    threadgroup float* bs, uint3 group, uint3 groups, uint tid) {
    constexpr int CN = BN / 64, BLOADS = BK * BN / 256;
    typedef LhsLoader<BM, 256, CONV> Loader;
    const int tx = tid % 16, ty = tid / 16;
    const int row0 = group.y * BM, col0 = group.x * BN;
    float acc[8][4 * CN];
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 4 * CN; j++) acc[i][j] = 0.0f;
    int kbegin, kend;
    split_range(p.k, BK, group, groups, kbegin, kend);
    Loader loader(a, p, row0, kbegin, tid);
    float av[Loader::ROWS], bv[BLOADS];

    // Fetch the first tile into registers, then threadgroup buffer 0.
    loader.fetch(av, p, kbegin, kend);
    for (int j = 0; j < BLOADS; j++) {
        int e = tid + 256 * j, kk = e / BN, col = e % BN;
        bv[j] = (kbegin + kk < kend && col0 + col < p.n) ? b[(size_t)(kbegin + kk) * p.n + col0 + col] : 0.0f;
    }
    for (int i = 0; i < Loader::ROWS; i++) as[(tid % 8) * AS + tid / 8 + 32 * i] = av[i];
    for (int j = 0; j < BLOADS; j++) {
        int e = tid + 256 * j;
        bs[(e / BN) * BN + e % BN] = bv[j];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    int buffer = 0;
    for (int k0 = kbegin; k0 < kend; k0 += BK) {
        const bool more = k0 + BK < kend;
        if (more) {
            loader.fetch(av, p, k0 + BK, kend);
            for (int j = 0; j < BLOADS; j++) {
                int e = tid + 256 * j, kk = e / BN, col = e % BN;
                bv[j] = (k0 + BK + kk < kend && col0 + col < p.n) ? b[(size_t)(k0 + BK + kk) * p.n + col0 + col] : 0.0f;
            }
        }
        threadgroup const float* at = as + buffer * BK * AS;
        threadgroup const float* bt = bs + buffer * BK * BN;
        for (int kk = 0; kk < BK; kk++) {
            float x[8], y[4 * CN];
            const float4 lo = *(threadgroup const float4*)(at + kk * AS + 4 * ty);
            const float4 hi = *(threadgroup const float4*)(at + kk * AS + 64 + 4 * ty);
            x[0] = lo.x; x[1] = lo.y; x[2] = lo.z; x[3] = lo.w;
            x[4] = hi.x; x[5] = hi.y; x[6] = hi.z; x[7] = hi.w;
            for (int j = 0; j < CN; j++) {
                const float4 c = *(threadgroup const float4*)(bt + kk * BN + 64 * j + 4 * tx);
                y[4 * j] = c.x; y[4 * j + 1] = c.y; y[4 * j + 2] = c.z; y[4 * j + 3] = c.w;
            }
            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 4 * CN; j++) acc[i][j] = fma(x[i], y[j], acc[i][j]);
        }
        if (more) {
            threadgroup float* an = as + (buffer ^ 1) * BK * AS;
            threadgroup float* bn = bs + (buffer ^ 1) * BK * BN;
            for (int i = 0; i < Loader::ROWS; i++) an[(tid % 8) * AS + tid / 8 + 32 * i] = av[i];
            for (int j = 0; j < BLOADS; j++) {
                int e = tid + 256 * j;
                bn[(e / BN) * BN + e % BN] = bv[j];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
            buffer ^= 1;
        }
    }
    for (int i = 0; i < 8; i++) {
        int row = row0 + (i < 4 ? 4 * ty + i : 64 + 4 * ty + i - 4);
        if (row >= p.m) continue;
        for (int j = 0; j < 4 * CN; j++) {
            int col = col0 + 64 * (j / 4) + 4 * tx + j % 4;
            if (col < p.n) store_sum(out, acc[i][j], row, col, p, bias, residual, group, groups);
        }
    }
}

#define PRODUCT(NAME, CALL, AS_LEN, BS_LEN)                                                               \
    kernel void NAME(device const float* a [[buffer(0)]], device const float* b [[buffer(1)]],          \
                     device float* out [[buffer(2)]], device const float* bias [[buffer(3)]],             \
                     device const float* residual [[buffer(4)]], constant ProductParams& p [[buffer(5)]], \
                     uint3 group [[threadgroup_position_in_grid]], uint3 groups [[threadgroups_per_grid]], \
                     uint tid [[thread_index_in_threadgroup]]) {                                          \
        threadgroup float as[AS_LEN], bs[BS_LEN];                                                        \
        CALL(a, b, out, bias, residual, p, as, bs, group, groups, tid);                                  \
    }

PRODUCT(gemm_128x64, (gemm128<64, false>), 2 * BK * AS, 2 * BK * 64)
PRODUCT(gemm_64x32, gemm_small<false>, 8 * 64, 8 * 32)
PRODUCT(conv_128x64, (gemm128<64, true>), 2 * BK * AS, 2 * BK * 64)

// out = act(sum of the `splits` partial sums in ws (each m x n) + bias), for split-K products.
struct ReduceParams {
    int m, n, ldc, splits, act;
};

kernel void splitk_reduce(device const float* ws [[buffer(0)]], device float* out [[buffer(1)]],
                          device const float* bias [[buffer(2)]], constant ReduceParams& p [[buffer(3)]],
                          uint e [[thread_position_in_grid]]) {
    if ((int)e >= p.m * p.n) return;
    int row = e / p.n, col = e % p.n;
    float v = 0.0f;
    for (int z = 0; z < p.splits; z++) v += ws[((size_t)z * p.m + row) * p.n + col];
    out[(size_t)row * p.ldc + col] = activate(v + bias[col], p.act);
}

inline float sigmoid(float x) { return 1.0f / (1.0f + precise::exp(-x)); }

// The cell update of step s of both directions (grid y) of a bidirectional LSTM layer with
// hidden size h over `steps` steps, after the gates' products: step s of the forward direction is
// position s, of the reverse one steps - 1 - s.
//   pre:  [n][steps][2][4][h] gate pre-activations (PyTorch's gate order i, f, g, o)
//   out:  [n][steps][2][h] hidden states
//   cell: [2][n][h]
struct CellParams {
    int n, steps, h, s;
};

kernel void lstm_cell(device const float* pre [[buffer(0)]], device float* out [[buffer(1)]],
                      device float* cell [[buffer(2)]], constant CellParams& p [[buffer(3)]],
                      uint2 at [[thread_position_in_grid]]) {
    const int h = p.h, d = at.y;
    if ((int)at.x >= p.n * h) return;
    const int row = at.x / h, u = at.x % h;
    const int t = d == 0 ? p.s : p.steps - 1 - p.s;
    device const float* g = pre + (((size_t)row * p.steps + t) * 2 + d) * 4 * h;
    device float* c = cell + ((size_t)d * p.n + row) * h + u;
    float gi = sigmoid(g[u]);
    float gf = sigmoid(g[h + u]);
    float gg = precise::tanh(g[2 * h + u]);
    float go = sigmoid(g[3 * h + u]);
    float cv = (p.s > 0 ? gf * *c : 0.0f) + gi * gg;
    *c = cv;
    out[((size_t)row * p.steps + t) * 2 * h + (size_t)d * h + u] = go * precise::tanh(cv);
}

// x[i] = (float)in[i] / divisor.
struct NormalizeParams {
    uint len;
    float divisor;
};

kernel void int8_to_float(device const char* in [[buffer(0)]], device float* x [[buffer(1)]],
                          constant NormalizeParams& p [[buffer(2)]], uint i [[thread_position_in_grid]]) {
    if (i < p.len) x[i] = (float)in[i] / p.divisor;
}

kernel void int32_to_float(device const int* in [[buffer(0)]], device float* x [[buffer(1)]],
                           constant uint& len [[buffer(2)]], uint i [[thread_position_in_grid]]) {
    if (i < len) x[i] = (float)in[i];
}

// PyramidPolling: TF-style "same" max pools of each [h][w][c] map to 3x3, 2x2 and 1x1 grids, the
// extra padding row or column at the end and counting as zero, each flattened [row][col][c].
// out: [images][14 c].
struct PoolParams {
    int images, h, w, c;
};

kernel void pyramid_pool(device const float* in [[buffer(0)]], device float* out [[buffer(1)]],
                         constant PoolParams& p [[buffer(2)]], uint e [[thread_position_in_grid]]) {
    const int c = p.c, h = p.h, w = p.w, features = 14 * c;
    if ((int)e >= p.images * features) return;
    int image = e / features, f = e % features;
    int size, cell;
    if (f < 9 * c) { size = 3; cell = f / c; }
    else if (f < 13 * c) { size = 2; cell = (f - 9 * c) / c; }
    else { size = 1; cell = 0; }
    int ch = f % c;
    int wh = (h + size - 1) / size, ww = (w + size - 1) / size;
    int out_h = (h + wh - 1) / wh, out_w = (w + ww - 1) / ww;
    int pad_h = max(out_h * wh - h, 0), pad_w = max(out_w * ww - w, 0);
    int oy = cell / out_w, ox = cell % out_w;
    int y0 = oy * wh - pad_h / 2, x0 = ox * ww - pad_w / 2;
    int y1 = y0 + wh, x1 = x0 + ww;
    bool padded = y0 < 0 || x0 < 0 || y1 > h || x1 > w;
    y0 = max(y0, 0); x0 = max(x0, 0);
    y1 = min(y1, h); x1 = min(x1, w);
    float m = padded ? 0.0f : -INFINITY;
    device const float* map = in + (size_t)image * h * w * c;
    for (int y = y0; y < y1; y++)
        for (int x = x0; x < x1; x++) m = max(m, map[((size_t)y * w + x) * c + ch]);
    out[(size_t)image * features + f] = m;
}

// Softmax over each row of logits [n][classes], written to out[row * ldo + offset + class].
struct SoftmaxParams {
    int n, classes, ldo, offset;
};

kernel void softmax_rows(device const float* logits [[buffer(0)]], device float* out [[buffer(1)]],
                         constant SoftmaxParams& p [[buffer(2)]], uint row [[thread_position_in_grid]]) {
    if ((int)row >= p.n) return;
    device const float* l = logits + (size_t)row * p.classes;
    float m = -INFINITY;
    for (int i = 0; i < p.classes; i++) m = max(m, l[i]);
    float sum = 0.0f;
    for (int i = 0; i < p.classes; i++) sum += precise::exp(l[i] - m);
    for (int i = 0; i < p.classes; i++) out[(size_t)row * p.ldo + p.offset + i] = precise::exp(l[i] - m) / sum;
}
