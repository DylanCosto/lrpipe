//! Port of the Striped Smith-Waterman library the original bundles (`preprocess/realign/ssw.c`
//! and `ssw_cpp.cpp`), as `StripedSmithWaterman::Aligner` with its default scores uses it.
//! The library is Mengyao Zhao's, copyright Boston College, under the MIT license: its notice is
//! in `lrpipe-snv/LICENSE-SSW`.
//!
//! The striped kernels are reproduced rather than replaced by a textbook dynamic programme:
//! their lazy-F loop does not update `E` after correcting `H`, so scores and end positions
//! depend on the SSE2 lane layout. Both the 8-bit kernel (tried first) and the 16-bit one (on
//! overflow) are reproduced, with their saturating arithmetic: on x86-64 with the same SSE2
//! instructions in the same order (`sse2`), on aarch64 with their NEON counterparts lane for lane
//! (`neon`), elsewhere emulated lane by lane (`lanes`), which the tests check them against.

/// `kBaseTranslation`: `A`/`a` 0, `C` 1, `G` 2, `T`/`U` 3, anything else 4.
fn translate(base: u8) -> i8 {
    match base {
        b'A' | b'a' | b'U' | b'u' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Alphabet size of the score matrix (`score_matrix_size_`).
const N: usize = 5;
/// `Aligner()` defaults: match 4, mismatch 6, gap open 8, gap extend 2.
const MATCH: i8 = 4;
const MISMATCH: i8 = 6;
const GAP_OPEN: u8 = 8;
const GAP_EXTEND: u8 = 2;

/// `BuildSwScoreMatrix`: `N` mismatches everything, itself included.
fn score_matrix() -> [i8; N * N] {
    let mut m = [-MISMATCH; N * N];
    for i in 0..4 {
        m[i * N + i] = MATCH;
    }
    m
}

/// The result `Aligner::Align` fills in (the fields the original reads).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Alignment {
    /// `sw_score`.
    pub score: u16,
    /// `ref_begin`, 0-based.
    pub ref_begin: i32,
    /// `ref_end`, 0-based.
    pub ref_end: i32,
    /// `query_begin`, 0-based.
    pub query_begin: i32,
    /// `query_end`, 0-based.
    pub query_end: i32,
    /// `cigar_string`, with `=` and `X` instead of `M`.
    pub cigar: String,
}

/// An `Aligner` with a reference set (`SetReferenceSequence`).
#[derive(Debug, Clone)]
pub struct Aligner {
    reference: Vec<i8>,
    matrix: [i8; N * N],
}

impl Aligner {
    /// `Aligner()` then `SetReferenceSequence(reference)`.
    pub fn new(reference: &[u8]) -> Self {
        Self { reference: reference.iter().map(|&b| translate(b)).collect(), matrix: score_matrix() }
    }

    /// `Aligner::Align(query, Filter(), &alignment)`, or `None` where it returns false (an empty
    /// query or reference).
    pub fn align(&self, query: &[u8]) -> Option<Alignment> {
        if self.reference.is_empty() || query.is_empty() {
            return None;
        }
        let query: Vec<i8> = query.iter().map(|&b| translate(b)).collect();
        let raw = ssw_align(&query, &self.reference, &self.matrix, query.len() as i32);
        let mut alignment = Alignment {
            score: raw.score1,
            ref_begin: raw.ref_begin1,
            ref_end: raw.ref_end1,
            query_begin: raw.read_begin1,
            query_end: raw.read_end1,
            cigar: String::new(),
        };
        alignment.cigar = match_mismatch_cigar(&alignment, &raw.cigar, &self.reference, &query);
        Some(alignment)
    }
}

/// `s_align`, without the second-best fields the original never reads.
struct RawAlignment {
    score1: u16,
    ref_begin1: i32,
    ref_end1: i32,
    read_begin1: i32,
    read_end1: i32,
    /// `(op, length)` of the banded traceback, `op` one of `M`, `I`, `D`.
    cigar: Vec<(u8, u32)>,
}

/// `alignment_end` (the best one).
struct End {
    score: u16,
    reference: i32,
    read: i32,
}

/// The query profile of `qP_byte` (16 lanes, biased) or `qP_word` (8 lanes).
fn profile<T: Copy>(read: &[i8], matrix: &[i8; N * N], lanes: usize, value: impl Fn(Option<i8>) -> T) -> Vec<T> {
    let read_len = read.len();
    let seg_len = read_len.div_ceil(lanes);
    let mut out = Vec::with_capacity(N * seg_len * lanes);
    for nt in 0..N {
        for i in 0..seg_len {
            for seg in 0..lanes {
                let j = i + seg * seg_len;
                out.push(value(read.get(j).map(|&r| matrix[nt * N + r as usize])));
            }
        }
    }
    out
}

/// The kernels emulated lane by lane: the reference for the SSE2 and NEON ones, and the kernels
/// on targets with neither.
#[cfg(any(test, not(any(target_arch = "x86_64", target_arch = "aarch64"))))]
mod lanes {
    use super::{End, GAP_EXTEND, GAP_OPEN};

    /// Lane-wise vectors of the 8-bit kernel.
    type V8 = [u8; 16];
    /// Lane-wise vectors of the 16-bit kernel.
    type V16 = [i16; 8];

    fn adds_u8(a: V8, b: V8) -> V8 {
        std::array::from_fn(|k| a[k].saturating_add(b[k]))
    }
    fn subs_u8(a: V8, b: V8) -> V8 {
        std::array::from_fn(|k| a[k].saturating_sub(b[k]))
    }
    fn max_u8(a: V8, b: V8) -> V8 {
        std::array::from_fn(|k| a[k].max(b[k]))
    }
    /// `_mm_slli_si128(v, 1)`: every lane moves up one, lane 0 becomes 0.
    fn shift_u8(a: V8) -> V8 {
        std::array::from_fn(|k| if k == 0 { 0 } else { a[k - 1] })
    }

    fn adds_i16(a: V16, b: V16) -> V16 {
        std::array::from_fn(|k| a[k].saturating_add(b[k]))
    }
    /// `_mm_subs_epu16`: unsigned saturating subtraction on the 16-bit patterns.
    fn subs_u16(a: V16, b: V16) -> V16 {
        std::array::from_fn(|k| (a[k] as u16).saturating_sub(b[k] as u16) as i16)
    }
    fn max_i16(a: V16, b: V16) -> V16 {
        std::array::from_fn(|k| a[k].max(b[k]))
    }
    /// `_mm_slli_si128(v, 2)`.
    fn shift_i16(a: V16) -> V16 {
        std::array::from_fn(|k| if k == 0 { 0 } else { a[k - 1] })
    }

    /// `sw_sse2_byte`: forward (`reverse = false`) or backward over `reference[..ref_len]`.
    pub(super) fn sw_byte(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[u8],
        terminate: u8,
        bias: u8,
    ) -> End {
        let seg_len = ((read_len + 15) / 16) as usize;
        let mut max: u8 = 0;
        let mut end_read = read_len - 1;
        let mut end_ref: i32 = -1;
        let mut store: Vec<V8> = vec![[0; 16]; seg_len];
        let mut load: Vec<V8> = vec![[0; 16]; seg_len];
        let mut pv_e: Vec<V8> = vec![[0; 16]; seg_len];
        let mut hmax: Vec<V8> = vec![[0; 16]; seg_len];
        let gap_o = [GAP_OPEN; 16];
        let gap_e = [GAP_EXTEND; 16];
        let v_bias = [bias; 16];
        let mut max_score: V8 = [0; 16];
        let mut max_mark: V8 = [0; 16];
        let (begin, end, step) = if reverse { (ref_len - 1, -1, -1) } else { (0, ref_len, 1) };
        let mut i = begin;
        while i != end {
            let mut v_f: V8 = [0; 16];
            let mut max_column: V8 = [0; 16];
            let mut v_h = shift_u8(store[seg_len - 1]);
            let base = reference[i as usize] as usize;
            let v_p = &profile[base * seg_len * 16..];
            std::mem::swap(&mut load, &mut store);
            for j in 0..seg_len {
                let p: V8 = v_p[j * 16..j * 16 + 16].try_into().expect("16 lanes");
                v_h = adds_u8(v_h, p);
                v_h = subs_u8(v_h, v_bias);
                let mut e = pv_e[j];
                v_h = max_u8(v_h, e);
                v_h = max_u8(v_h, v_f);
                max_column = max_u8(max_column, v_h);
                store[j] = v_h;
                v_h = subs_u8(v_h, gap_o);
                e = subs_u8(e, gap_e);
                e = max_u8(e, v_h);
                pv_e[j] = e;
                v_f = subs_u8(v_f, gap_e);
                v_f = max_u8(v_f, v_h);
                v_h = load[j];
            }
            // Lazy F loop.
            let mut j = 0;
            let mut v_h = store[j];
            v_f = shift_u8(v_f);
            let all_done = |v_h: V8, v_f: V8| subs_u8(v_f, subs_u8(v_h, gap_o)).iter().all(|&x| x == 0);
            while !all_done(v_h, v_f) {
                v_h = max_u8(v_h, v_f);
                max_column = max_u8(max_column, v_h);
                store[j] = v_h;
                v_f = subs_u8(v_f, gap_e);
                j += 1;
                if j >= seg_len {
                    j = 0;
                    v_f = shift_u8(v_f);
                }
                v_h = store[j];
            }
            max_score = max_u8(max_score, max_column);
            if max_mark != max_score {
                max_mark = max_score;
                let temp = *max_score.iter().max().expect("16 lanes");
                if temp > max {
                    max = temp;
                    if u32::from(max) + u32::from(bias) >= 255 {
                        break;
                    }
                    end_ref = i;
                    hmax.clone_from(&store);
                }
            }
            let column_max = *max_column.iter().max().expect("16 lanes");
            if column_max == terminate {
                break;
            }
            i += step;
        }
        for (index, &value) in hmax.iter().flatten().enumerate() {
            if value == max {
                let temp = (index / 16 + index % 16 * seg_len) as i32;
                end_read = end_read.min(temp);
            }
        }
        let score = if u32::from(max) + u32::from(bias) >= 255 { 255 } else { u16::from(max) };
        End { score, reference: end_ref, read: end_read }
    }

    /// `sw_sse2_word`.
    pub(super) fn sw_word(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[i16],
        terminate: u16,
    ) -> End {
        let seg_len = ((read_len + 7) / 8) as usize;
        let mut max: u16 = 0;
        let mut end_read = read_len - 1;
        let mut end_ref: i32 = 0;
        let mut store: Vec<V16> = vec![[0; 8]; seg_len];
        let mut load: Vec<V16> = vec![[0; 8]; seg_len];
        let mut pv_e: Vec<V16> = vec![[0; 8]; seg_len];
        let mut hmax: Vec<V16> = vec![[0; 8]; seg_len];
        let gap_o = [i16::from(GAP_OPEN); 8];
        let gap_e = [i16::from(GAP_EXTEND); 8];
        let mut max_score: V16 = [0; 8];
        let mut max_mark: V16 = [0; 8];
        let (begin, end, step) = if reverse { (ref_len - 1, -1, -1) } else { (0, ref_len, 1) };
        let mut i = begin;
        while i != end {
            let mut v_f: V16 = [0; 8];
            let mut v_h = shift_i16(store[seg_len - 1]);
            let mut max_column: V16 = [0; 8];
            let base = reference[i as usize] as usize;
            let v_p = &profile[base * seg_len * 8..];
            std::mem::swap(&mut load, &mut store);
            for j in 0..seg_len {
                let p: V16 = v_p[j * 8..j * 8 + 8].try_into().expect("8 lanes");
                v_h = adds_i16(v_h, p);
                let mut e = pv_e[j];
                v_h = max_i16(v_h, e);
                v_h = max_i16(v_h, v_f);
                max_column = max_i16(max_column, v_h);
                store[j] = v_h;
                v_h = subs_u16(v_h, gap_o);
                e = subs_u16(e, gap_e);
                e = max_i16(e, v_h);
                pv_e[j] = e;
                v_f = subs_u16(v_f, gap_e);
                v_f = max_i16(v_f, v_h);
                v_h = load[j];
            }
            // Lazy F loop.
            'lazy: for _ in 0..8 {
                v_f = shift_i16(v_f);
                for slot in store.iter_mut() {
                    let mut v_h = max_i16(*slot, v_f);
                    max_column = max_i16(max_column, v_h);
                    *slot = v_h;
                    v_h = subs_u16(v_h, gap_o);
                    v_f = subs_u16(v_f, gap_e);
                    if !v_f.iter().zip(&v_h).any(|(f, h)| f > h) {
                        break 'lazy;
                    }
                }
            }
            max_score = max_i16(max_score, max_column);
            if max_mark != max_score {
                max_mark = max_score;
                let temp = *max_score.iter().max().expect("8 lanes") as u16;
                if temp > max {
                    max = temp;
                    end_ref = i;
                    hmax.clone_from(&store);
                }
            }
            let column_max = *max_column.iter().max().expect("8 lanes") as u16;
            if column_max == terminate {
                break;
            }
            i += step;
        }
        for (index, &value) in hmax.iter().flatten().enumerate() {
            if value as u16 == max {
                let temp = (index / 8 + index % 8 * seg_len) as i32;
                end_read = end_read.min(temp);
            }
        }
        End { score: max, reference: end_ref, read: end_read }
    }
}

/// `sw_sse2_byte` and `sw_sse2_word` with SSE2 intrinsics, instruction for instruction (x86-64
/// always has SSE2).
#[cfg(target_arch = "x86_64")]
mod sse2 {
    use super::{End, GAP_EXTEND, GAP_OPEN};
    use std::arch::x86_64::{
        __m128i, _mm_adds_epi16, _mm_adds_epu8, _mm_cmpeq_epi8, _mm_cmpeq_epi16, _mm_cmpgt_epi16, _mm_cvtsi128_si32,
        _mm_loadu_si128, _mm_max_epi16, _mm_max_epu8, _mm_movemask_epi8, _mm_set1_epi8, _mm_set1_epi16,
        _mm_setzero_si128, _mm_slli_si128, _mm_srli_si128, _mm_storeu_si128, _mm_subs_epu8, _mm_subs_epu16,
    };

    /// The query profile as vectors, `lanes` elements each.
    #[target_feature(enable = "sse2")]
    fn vectors<T>(profile: &[T], lanes: usize) -> Vec<__m128i> {
        assert_eq!(size_of::<T>() * lanes, 16);
        profile
            .chunks_exact(lanes)
            // SAFETY: each chunk holds 16 bytes; the load is unaligned.
            .map(|c| unsafe { _mm_loadu_si128(c.as_ptr().cast::<__m128i>()) })
            .collect()
    }

    #[target_feature(enable = "sse2")]
    fn bytes(v: __m128i) -> [u8; 16] {
        let mut out = [0u8; 16];
        // SAFETY: `out` holds 16 bytes; the store is unaligned.
        unsafe { _mm_storeu_si128(out.as_mut_ptr().cast::<__m128i>(), v) };
        out
    }

    #[target_feature(enable = "sse2")]
    fn words(v: __m128i) -> [i16; 8] {
        let mut out = [0i16; 8];
        // SAFETY: `out` holds 16 bytes; the store is unaligned.
        unsafe { _mm_storeu_si128(out.as_mut_ptr().cast::<__m128i>(), v) };
        out
    }

    /// `max16`: the largest unsigned byte.
    #[target_feature(enable = "sse2")]
    fn max_u8(v: __m128i) -> u8 {
        let v = _mm_max_epu8(v, _mm_srli_si128::<8>(v));
        let v = _mm_max_epu8(v, _mm_srli_si128::<4>(v));
        let v = _mm_max_epu8(v, _mm_srli_si128::<2>(v));
        let v = _mm_max_epu8(v, _mm_srli_si128::<1>(v));
        _mm_cvtsi128_si32(v) as u8
    }

    /// `max8`: the largest signed word.
    #[target_feature(enable = "sse2")]
    fn max_i16(v: __m128i) -> i16 {
        let v = _mm_max_epi16(v, _mm_srli_si128::<8>(v));
        let v = _mm_max_epi16(v, _mm_srli_si128::<4>(v));
        let v = _mm_max_epi16(v, _mm_srli_si128::<2>(v));
        _mm_cvtsi128_si32(v) as i16
    }

    /// `sw_sse2_byte`: forward (`reverse = false`) or backward over `reference[..ref_len]`.
    pub(super) fn sw_byte(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[u8],
        terminate: u8,
        bias: u8,
    ) -> End {
        // SAFETY: x86-64 always has SSE2.
        unsafe { byte_kernel(reference, reverse, ref_len, read_len, profile, terminate, bias) }
    }

    #[target_feature(enable = "sse2")]
    fn byte_kernel(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[u8],
        terminate: u8,
        bias: u8,
    ) -> End {
        let seg_len = ((read_len + 15) / 16) as usize;
        let profile = vectors(profile, 16);
        let zero = _mm_setzero_si128();
        let mut max: u8 = 0;
        let mut end_read = read_len - 1;
        let mut end_ref: i32 = -1;
        let mut store = vec![zero; seg_len];
        let mut load = vec![zero; seg_len];
        let mut pv_e = vec![zero; seg_len];
        let mut hmax = vec![zero; seg_len];
        let gap_o = _mm_set1_epi8(GAP_OPEN as i8);
        let gap_e = _mm_set1_epi8(GAP_EXTEND as i8);
        let v_bias = _mm_set1_epi8(bias as i8);
        let mut max_score = zero;
        let mut max_mark = zero;
        let (begin, end, step) = if reverse { (ref_len - 1, -1, -1) } else { (0, ref_len, 1) };
        let mut i = begin;
        while i != end {
            let mut v_f = zero;
            let mut max_column = zero;
            let mut v_h = _mm_slli_si128::<1>(store[seg_len - 1]);
            let base = reference[i as usize] as usize;
            let v_p = &profile[base * seg_len..(base + 1) * seg_len];
            std::mem::swap(&mut load, &mut store);
            for j in 0..seg_len {
                v_h = _mm_adds_epu8(v_h, v_p[j]);
                v_h = _mm_subs_epu8(v_h, v_bias);
                let mut e = pv_e[j];
                v_h = _mm_max_epu8(v_h, e);
                v_h = _mm_max_epu8(v_h, v_f);
                max_column = _mm_max_epu8(max_column, v_h);
                store[j] = v_h;
                v_h = _mm_subs_epu8(v_h, gap_o);
                e = _mm_subs_epu8(e, gap_e);
                e = _mm_max_epu8(e, v_h);
                pv_e[j] = e;
                v_f = _mm_subs_epu8(v_f, gap_e);
                v_f = _mm_max_epu8(v_f, v_h);
                v_h = load[j];
            }
            // Lazy F loop.
            let mut j = 0;
            let mut v_h = store[j];
            v_f = _mm_slli_si128::<1>(v_f);
            let all_done = |v_h: __m128i, v_f: __m128i| {
                _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(v_f, _mm_subs_epu8(v_h, gap_o)), zero)) == 0xffff
            };
            while !all_done(v_h, v_f) {
                v_h = _mm_max_epu8(v_h, v_f);
                max_column = _mm_max_epu8(max_column, v_h);
                store[j] = v_h;
                v_f = _mm_subs_epu8(v_f, gap_e);
                j += 1;
                if j >= seg_len {
                    j = 0;
                    v_f = _mm_slli_si128::<1>(v_f);
                }
                v_h = store[j];
            }
            max_score = _mm_max_epu8(max_score, max_column);
            if _mm_movemask_epi8(_mm_cmpeq_epi8(max_mark, max_score)) != 0xffff {
                max_mark = max_score;
                let temp = max_u8(max_score);
                if temp > max {
                    max = temp;
                    if u32::from(max) + u32::from(bias) >= 255 {
                        break;
                    }
                    end_ref = i;
                    hmax.copy_from_slice(&store);
                }
            }
            if max_u8(max_column) == terminate {
                break;
            }
            i += step;
        }
        for (index, value) in hmax.iter().flat_map(|&v| bytes(v)).enumerate() {
            if value == max {
                let temp = (index / 16 + index % 16 * seg_len) as i32;
                end_read = end_read.min(temp);
            }
        }
        let score = if u32::from(max) + u32::from(bias) >= 255 { 255 } else { u16::from(max) };
        End { score, reference: end_ref, read: end_read }
    }

    /// `sw_sse2_word`.
    pub(super) fn sw_word(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[i16],
        terminate: u16,
    ) -> End {
        // SAFETY: x86-64 always has SSE2.
        unsafe { word_kernel(reference, reverse, ref_len, read_len, profile, terminate) }
    }

    #[target_feature(enable = "sse2")]
    fn word_kernel(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[i16],
        terminate: u16,
    ) -> End {
        let seg_len = ((read_len + 7) / 8) as usize;
        let profile = vectors(profile, 8);
        let zero = _mm_setzero_si128();
        let mut max: u16 = 0;
        let mut end_read = read_len - 1;
        let mut end_ref: i32 = 0;
        let mut store = vec![zero; seg_len];
        let mut load = vec![zero; seg_len];
        let mut pv_e = vec![zero; seg_len];
        let mut hmax = vec![zero; seg_len];
        let gap_o = _mm_set1_epi16(i16::from(GAP_OPEN));
        let gap_e = _mm_set1_epi16(i16::from(GAP_EXTEND));
        let mut max_score = zero;
        let mut max_mark = zero;
        let (begin, end, step) = if reverse { (ref_len - 1, -1, -1) } else { (0, ref_len, 1) };
        let mut i = begin;
        while i != end {
            let mut v_f = zero;
            let mut v_h = _mm_slli_si128::<2>(store[seg_len - 1]);
            let mut max_column = zero;
            let base = reference[i as usize] as usize;
            let v_p = &profile[base * seg_len..(base + 1) * seg_len];
            std::mem::swap(&mut load, &mut store);
            for j in 0..seg_len {
                v_h = _mm_adds_epi16(v_h, v_p[j]);
                let mut e = pv_e[j];
                v_h = _mm_max_epi16(v_h, e);
                v_h = _mm_max_epi16(v_h, v_f);
                max_column = _mm_max_epi16(max_column, v_h);
                store[j] = v_h;
                v_h = _mm_subs_epu16(v_h, gap_o);
                e = _mm_subs_epu16(e, gap_e);
                e = _mm_max_epi16(e, v_h);
                pv_e[j] = e;
                v_f = _mm_subs_epu16(v_f, gap_e);
                v_f = _mm_max_epi16(v_f, v_h);
                v_h = load[j];
            }
            // Lazy F loop.
            'lazy: for _ in 0..8 {
                v_f = _mm_slli_si128::<2>(v_f);
                for slot in store.iter_mut() {
                    let mut v_h = _mm_max_epi16(*slot, v_f);
                    max_column = _mm_max_epi16(max_column, v_h);
                    *slot = v_h;
                    v_h = _mm_subs_epu16(v_h, gap_o);
                    v_f = _mm_subs_epu16(v_f, gap_e);
                    if _mm_movemask_epi8(_mm_cmpgt_epi16(v_f, v_h)) == 0 {
                        break 'lazy;
                    }
                }
            }
            max_score = _mm_max_epi16(max_score, max_column);
            if _mm_movemask_epi8(_mm_cmpeq_epi16(max_mark, max_score)) != 0xffff {
                max_mark = max_score;
                let temp = max_i16(max_score) as u16;
                if temp > max {
                    max = temp;
                    end_ref = i;
                    hmax.copy_from_slice(&store);
                }
            }
            if max_i16(max_column) as u16 == terminate {
                break;
            }
            i += step;
        }
        for (index, value) in hmax.iter().flat_map(|&v| words(v)).enumerate() {
            if value as u16 == max {
                let temp = (index / 8 + index % 8 * seg_len) as i32;
                end_read = end_read.min(temp);
            }
        }
        End { score: max, reference: end_ref, read: end_read }
    }
}

/// `sw_sse2_byte` and `sw_sse2_word` with NEON (aarch64 always has it): each SSE2 instruction's
/// counterpart on the same 16 lanes, so the same values in the same order. `_mm_slli_si128` is
/// `vext` with a zero vector, the `movemask` tests are across-vector maxima and minima.
#[cfg(target_arch = "aarch64")]
mod neon {
    use super::{End, GAP_EXTEND, GAP_OPEN};
    use std::arch::aarch64::{
        int16x8_t, uint8x16_t, vceqq_s16, vceqq_u8, vcgtq_s16, vdupq_n_s16, vdupq_n_u8, vextq_s16, vextq_u8, vld1q_s16,
        vld1q_u8, vmaxq_s16, vmaxq_u8, vmaxvq_s16, vmaxvq_u8, vmaxvq_u16, vminvq_u8, vminvq_u16, vqaddq_s16, vqaddq_u8,
        vqsubq_u8, vqsubq_u16, vreinterpretq_s16_u16, vreinterpretq_u16_s16, vst1q_s16, vst1q_u8,
    };

    /// `_mm_slli_si128(v, 1)`: every lane moves up one, lane 0 becomes 0.
    #[target_feature(enable = "neon")]
    fn shift_u8(v: uint8x16_t) -> uint8x16_t {
        vextq_u8::<15>(vdupq_n_u8(0), v)
    }

    /// `_mm_slli_si128(v, 2)`.
    #[target_feature(enable = "neon")]
    fn shift_i16(v: int16x8_t) -> int16x8_t {
        vextq_s16::<7>(vdupq_n_s16(0), v)
    }

    /// `_mm_subs_epu16`: unsigned saturating subtraction on the 16-bit patterns.
    #[target_feature(enable = "neon")]
    fn subs_u16(a: int16x8_t, b: int16x8_t) -> int16x8_t {
        vreinterpretq_s16_u16(vqsubq_u16(vreinterpretq_u16_s16(a), vreinterpretq_u16_s16(b)))
    }

    /// Whether every lane of `a` equals `b`'s.
    #[target_feature(enable = "neon")]
    fn same_u8(a: uint8x16_t, b: uint8x16_t) -> bool {
        vminvq_u8(vceqq_u8(a, b)) == u8::MAX
    }

    #[target_feature(enable = "neon")]
    fn same_i16(a: int16x8_t, b: int16x8_t) -> bool {
        vminvq_u16(vceqq_s16(a, b)) == u16::MAX
    }

    #[target_feature(enable = "neon")]
    fn bytes(v: uint8x16_t) -> [u8; 16] {
        let mut out = [0u8; 16];
        // SAFETY: `out` holds 16 bytes.
        unsafe { vst1q_u8(out.as_mut_ptr(), v) };
        out
    }

    #[target_feature(enable = "neon")]
    fn words(v: int16x8_t) -> [i16; 8] {
        let mut out = [0i16; 8];
        // SAFETY: `out` holds 8 words.
        unsafe { vst1q_s16(out.as_mut_ptr(), v) };
        out
    }

    /// `sw_sse2_byte`: forward (`reverse = false`) or backward over `reference[..ref_len]`.
    pub(super) fn sw_byte(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[u8],
        terminate: u8,
        bias: u8,
    ) -> End {
        // SAFETY: aarch64 always has NEON.
        unsafe { byte_kernel(reference, reverse, ref_len, read_len, profile, terminate, bias) }
    }

    #[target_feature(enable = "neon")]
    fn byte_kernel(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[u8],
        terminate: u8,
        bias: u8,
    ) -> End {
        let seg_len = ((read_len + 15) / 16) as usize;
        // SAFETY: each chunk holds 16 bytes.
        let profile: Vec<uint8x16_t> =
            profile.as_chunks::<16>().0.iter().map(|c| unsafe { vld1q_u8(c.as_ptr()) }).collect();
        let zero = vdupq_n_u8(0);
        let mut max: u8 = 0;
        let mut end_read = read_len - 1;
        let mut end_ref: i32 = -1;
        let mut store = vec![zero; seg_len];
        let mut load = vec![zero; seg_len];
        let mut pv_e = vec![zero; seg_len];
        let mut hmax = vec![zero; seg_len];
        let gap_o = vdupq_n_u8(GAP_OPEN);
        let gap_e = vdupq_n_u8(GAP_EXTEND);
        let v_bias = vdupq_n_u8(bias);
        let mut max_score = zero;
        let mut max_mark = zero;
        let (begin, end, step) = if reverse { (ref_len - 1, -1, -1) } else { (0, ref_len, 1) };
        let mut i = begin;
        while i != end {
            let mut v_f = zero;
            let mut max_column = zero;
            let mut v_h = shift_u8(store[seg_len - 1]);
            let base = reference[i as usize] as usize;
            let v_p = &profile[base * seg_len..(base + 1) * seg_len];
            std::mem::swap(&mut load, &mut store);
            for j in 0..seg_len {
                v_h = vqaddq_u8(v_h, v_p[j]);
                v_h = vqsubq_u8(v_h, v_bias);
                let mut e = pv_e[j];
                v_h = vmaxq_u8(v_h, e);
                v_h = vmaxq_u8(v_h, v_f);
                max_column = vmaxq_u8(max_column, v_h);
                store[j] = v_h;
                v_h = vqsubq_u8(v_h, gap_o);
                e = vqsubq_u8(e, gap_e);
                e = vmaxq_u8(e, v_h);
                pv_e[j] = e;
                v_f = vqsubq_u8(v_f, gap_e);
                v_f = vmaxq_u8(v_f, v_h);
                v_h = load[j];
            }
            // Lazy F loop.
            let mut j = 0;
            let mut v_h = store[j];
            v_f = shift_u8(v_f);
            let all_done = |v_h: uint8x16_t, v_f: uint8x16_t| vmaxvq_u8(vqsubq_u8(v_f, vqsubq_u8(v_h, gap_o))) == 0;
            while !all_done(v_h, v_f) {
                v_h = vmaxq_u8(v_h, v_f);
                max_column = vmaxq_u8(max_column, v_h);
                store[j] = v_h;
                v_f = vqsubq_u8(v_f, gap_e);
                j += 1;
                if j >= seg_len {
                    j = 0;
                    v_f = shift_u8(v_f);
                }
                v_h = store[j];
            }
            max_score = vmaxq_u8(max_score, max_column);
            if !same_u8(max_mark, max_score) {
                max_mark = max_score;
                let temp = vmaxvq_u8(max_score);
                if temp > max {
                    max = temp;
                    if u32::from(max) + u32::from(bias) >= 255 {
                        break;
                    }
                    end_ref = i;
                    hmax.copy_from_slice(&store);
                }
            }
            if vmaxvq_u8(max_column) == terminate {
                break;
            }
            i += step;
        }
        for (index, value) in hmax.iter().flat_map(|&v| bytes(v)).enumerate() {
            if value == max {
                let temp = (index / 16 + index % 16 * seg_len) as i32;
                end_read = end_read.min(temp);
            }
        }
        let score = if u32::from(max) + u32::from(bias) >= 255 { 255 } else { u16::from(max) };
        End { score, reference: end_ref, read: end_read }
    }

    /// `sw_sse2_word`.
    pub(super) fn sw_word(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[i16],
        terminate: u16,
    ) -> End {
        // SAFETY: aarch64 always has NEON.
        unsafe { word_kernel(reference, reverse, ref_len, read_len, profile, terminate) }
    }

    #[target_feature(enable = "neon")]
    fn word_kernel(
        reference: &[i8],
        reverse: bool,
        ref_len: i32,
        read_len: i32,
        profile: &[i16],
        terminate: u16,
    ) -> End {
        let seg_len = ((read_len + 7) / 8) as usize;
        // SAFETY: each chunk holds 8 words.
        let profile: Vec<int16x8_t> =
            profile.as_chunks::<8>().0.iter().map(|c| unsafe { vld1q_s16(c.as_ptr()) }).collect();
        let zero = vdupq_n_s16(0);
        let mut max: u16 = 0;
        let mut end_read = read_len - 1;
        let mut end_ref: i32 = 0;
        let mut store = vec![zero; seg_len];
        let mut load = vec![zero; seg_len];
        let mut pv_e = vec![zero; seg_len];
        let mut hmax = vec![zero; seg_len];
        let gap_o = vdupq_n_s16(i16::from(GAP_OPEN));
        let gap_e = vdupq_n_s16(i16::from(GAP_EXTEND));
        let mut max_score = zero;
        let mut max_mark = zero;
        let (begin, end, step) = if reverse { (ref_len - 1, -1, -1) } else { (0, ref_len, 1) };
        let mut i = begin;
        while i != end {
            let mut v_f = zero;
            let mut v_h = shift_i16(store[seg_len - 1]);
            let mut max_column = zero;
            let base = reference[i as usize] as usize;
            let v_p = &profile[base * seg_len..(base + 1) * seg_len];
            std::mem::swap(&mut load, &mut store);
            for j in 0..seg_len {
                v_h = vqaddq_s16(v_h, v_p[j]);
                let mut e = pv_e[j];
                v_h = vmaxq_s16(v_h, e);
                v_h = vmaxq_s16(v_h, v_f);
                max_column = vmaxq_s16(max_column, v_h);
                store[j] = v_h;
                v_h = subs_u16(v_h, gap_o);
                e = subs_u16(e, gap_e);
                e = vmaxq_s16(e, v_h);
                pv_e[j] = e;
                v_f = subs_u16(v_f, gap_e);
                v_f = vmaxq_s16(v_f, v_h);
                v_h = load[j];
            }
            // Lazy F loop.
            'lazy: for _ in 0..8 {
                v_f = shift_i16(v_f);
                for slot in store.iter_mut() {
                    let mut v_h = vmaxq_s16(*slot, v_f);
                    max_column = vmaxq_s16(max_column, v_h);
                    *slot = v_h;
                    v_h = subs_u16(v_h, gap_o);
                    v_f = subs_u16(v_f, gap_e);
                    if vmaxvq_u16(vcgtq_s16(v_f, v_h)) == 0 {
                        break 'lazy;
                    }
                }
            }
            max_score = vmaxq_s16(max_score, max_column);
            if !same_i16(max_mark, max_score) {
                max_mark = max_score;
                let temp = vmaxvq_s16(max_score) as u16;
                if temp > max {
                    max = temp;
                    end_ref = i;
                    hmax.copy_from_slice(&store);
                }
            }
            if vmaxvq_s16(max_column) as u16 == terminate {
                break;
            }
            i += step;
        }
        for (index, value) in hmax.iter().flat_map(|&v| words(v)).enumerate() {
            if value as u16 == max {
                let temp = (index / 8 + index % 8 * seg_len) as i32;
                end_read = end_read.min(temp);
            }
        }
        End { score: max, reference: end_ref, read: end_read }
    }
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
use lanes::{sw_byte, sw_word};
#[cfg(target_arch = "aarch64")]
use neon::{sw_byte, sw_word};
#[cfg(target_arch = "x86_64")]
use sse2::{sw_byte, sw_word};

/// `ssw_init(read, ..., score_size = 2)` then `ssw_align(..., flag = 0x0f, 0, 32767, maskLen)`.
fn ssw_align(read: &[i8], reference: &[i8], matrix: &[i8; N * N], _mask_len: i32) -> RawAlignment {
    let bias = matrix.iter().map(|&v| v.min(0)).min().unwrap_or(0).unsigned_abs();
    let ref_len = reference.len() as i32;
    let read_len = read.len() as i32;
    let byte_profile =
        |r: &[i8]| profile(r, matrix, 16, |v| v.map_or(bias, |s| (i16::from(s) + i16::from(bias)) as u8));
    let word_profile = |r: &[i8]| profile(r, matrix, 8, |v| v.map_or(0, i16::from));
    let mut word = false;
    let mut best = sw_byte(reference, false, ref_len, read_len, &byte_profile(read), u8::MAX, bias);
    if best.score == 255 {
        best = sw_word(reference, false, ref_len, read_len, &word_profile(read), u16::MAX);
        word = true;
    }
    let (score1, ref_end1, read_end1) = (best.score, best.reference, best.read);

    // The beginning: align the reversed read prefix backwards from the end.
    let reversed: Vec<i8> = (0..=read_end1).rev().map(|i| read[i as usize]).collect();
    let back = if word {
        sw_word(reference, true, ref_end1 + 1, read_end1 + 1, &word_profile(&reversed), score1)
    } else {
        sw_byte(reference, true, ref_end1 + 1, read_end1 + 1, &byte_profile(&reversed), score1 as u8, bias)
    };
    let ref_begin1 = back.reference;
    let read_begin1 = read_end1 - back.read;

    if ref_begin1 < 0 || read_begin1 < 0 {
        // Nothing aligned (score 0): the C code runs the banded alignment on `ref - 1`, one
        // base before the reference, which gives a single `M`.
        return RawAlignment { score1, ref_begin1, ref_end1, read_begin1, read_end1, cigar: vec![(b'M', 1)] };
    }
    let sub_ref = (ref_end1 - ref_begin1 + 1) as usize;
    let sub_read = (read_end1 - read_begin1 + 1) as usize;
    let band_width = (sub_ref as i32 - sub_read as i32).abs() + 1;
    let cigar = banded_sw(
        &reference[ref_begin1 as usize..],
        &read[read_begin1 as usize..],
        sub_ref,
        sub_read,
        i32::from(score1),
        band_width,
        matrix,
    );
    RawAlignment { score1, ref_begin1, ref_end1, read_begin1, read_end1, cigar }
}

/// `banded_sw`: the traceback of the best alignment, widening the band until it reaches
/// `score`.
fn banded_sw(
    reference: &[i8],
    read: &[i8],
    ref_len: usize,
    read_len: usize,
    score: i32,
    mut band_width: i32,
    matrix: &[i8; N * N],
) -> Vec<(u8, u32)> {
    let (gap_o, gap_e) = (i32::from(GAP_OPEN), i32::from(GAP_EXTEND));
    let mut direction: Vec<i8>;
    let mut width_d;
    // Index of `(i, j)` in a band row (`set_u`), and in a direction row (`set_d`).
    let set_u = |w: i32, i: i32, j: i32| (j - (i - w).max(0) + 1) as usize;
    let set_d = |w: i32, i: i32, j: i32, p: i32| ((j - (i - w).max(0)) * 3 + p) as usize;
    loop {
        let width = (band_width * 2 + 3) as usize;
        width_d = (band_width * 2 + 1) as usize;
        let mut h_b = vec![0i32; width + 1];
        let mut e_b = vec![0i32; width + 1];
        let mut h_c = vec![0i32; width + 1];
        direction = vec![0i8; width_d * read_len * 3];
        let mut max = 0;
        for i in 0..read_len as i32 {
            let beg = (i - band_width).max(0);
            let end = (ref_len as i32 - 1).min(i + band_width);
            let edge = ((end + 1) as usize).min(width - 1);
            let mut f = 0;
            h_b[0] = 0;
            e_b[0] = 0;
            h_b[edge] = 0;
            e_b[edge] = 0;
            h_c[0] = 0;
            let line = width_d * i as usize * 3;
            let mut u = 0;
            for j in beg..=end {
                u = set_u(band_width, i, j);
                let e = set_u(band_width, i - 1, j);
                let b = set_u(band_width, i, j - 1);
                let d = set_u(band_width, i - 1, j - 1);
                let de = line + set_d(band_width, i, j, 0);
                let df = line + set_d(band_width, i, j, 1);
                let dh = line + set_d(band_width, i, j, 2);

                let temp1 = if i == 0 { -gap_o } else { h_b[e] - gap_o };
                let temp2 = if i == 0 { -gap_e } else { e_b[e] - gap_e };
                e_b[u] = temp1.max(temp2);
                direction[de] = if temp1 > temp2 { 3 } else { 2 };

                let temp1 = h_c[b] - gap_o;
                let temp2 = f - gap_e;
                f = temp1.max(temp2);
                direction[df] = if temp1 > temp2 { 5 } else { 4 };

                let e1 = e_b[u].max(0);
                let f1 = f.max(0);
                let temp1 = e1.max(f1);
                let temp2 = h_b[d] + i32::from(matrix[reference[j as usize] as usize * N + read[i as usize] as usize]);
                h_c[u] = temp1.max(temp2);
                max = max.max(h_c[u]);
                direction[dh] = if temp1 <= temp2 {
                    1
                } else if e1 > f1 {
                    direction[de]
                } else {
                    direction[df]
                };
            }
            h_b[1..=u].copy_from_slice(&h_c[1..=u]);
        }
        band_width *= 2;
        if max >= score {
            break;
        }
    }
    band_width /= 2;

    // Traceback from the end of both sequences.
    let (mut i, mut j) = (read_len as i32 - 1, ref_len as i32 - 1);
    let mut e = 0u32;
    let mut line = width_d * i as usize * 3;
    let (mut op, mut prev_op) = (b'M', b'M');
    let mut matrix_id = 2;
    let mut cigar: Vec<(u8, u32)> = Vec::new();
    while i > 0 {
        match direction[line + set_d(band_width, i, j, matrix_id)] {
            1 => {
                i -= 1;
                j -= 1;
                matrix_id = 2;
                line -= width_d * 3;
                op = b'M';
            }
            2 => {
                i -= 1;
                matrix_id = 0;
                line -= width_d * 3;
                op = b'I';
            }
            3 => {
                i -= 1;
                matrix_id = 2;
                line -= width_d * 3;
                op = b'I';
            }
            4 => {
                j -= 1;
                matrix_id = 1;
                op = b'D';
            }
            5 => {
                j -= 1;
                matrix_id = 2;
                op = b'D';
            }
            other => unreachable!("SSW trace back error: {other}"),
        }
        if op == prev_op {
            e += 1;
        } else {
            cigar.push((prev_op, e));
            prev_op = op;
            e = 1;
        }
    }
    if op == b'M' {
        cigar.push((op, e + 1));
    } else {
        cigar.push((op, e));
        cigar.push((b'M', 1));
    }
    cigar.reverse();
    cigar
}

/// `ConvertAlignment` then `CalculateNumberMismatch`: soft clips around the traceback, with
/// each `M` split into `=` and `X` runs by comparing the translated bases.
fn match_mismatch_cigar(alignment: &Alignment, ops: &[(u8, u32)], reference: &[i8], query: &[i8]) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    if ops.is_empty() {
        return out;
    }
    if alignment.query_begin > 0 {
        write!(out, "{}S", alignment.query_begin).expect("writing to a String cannot fail");
    }
    // Signed: `ref_begin` is -1 when nothing aligned, and the C code then compares a byte
    // before the reference, which is taken as a mismatch here.
    let mut r = i64::from(alignment.ref_begin);
    let mut q = i64::from(alignment.query_begin);
    let base = |seq: &[i8], i: i64| usize::try_from(i).ok().and_then(|i| seq.get(i)).copied();
    let (mut in_m, mut in_x) = (false, false);
    let (mut length_m, mut length_x) = (0u32, 0u32);
    let flush = |out: &mut String, in_m: &mut bool, in_x: &mut bool, length_m: &mut u32, length_x: &mut u32| {
        if *in_m {
            write!(out, "{length_m}=").expect("writing to a String cannot fail");
        } else if *in_x {
            write!(out, "{length_x}X").expect("writing to a String cannot fail");
        }
        (*in_m, *in_x, *length_m, *length_x) = (false, false, 0, 0);
    };
    for &(op, length) in ops {
        match op {
            b'M' => {
                for _ in 0..length {
                    let (rb, qb) = (base(reference, r), base(query, q));
                    if rb.is_none() || rb != qb {
                        if in_m {
                            write!(out, "{length_m}=").expect("writing to a String cannot fail");
                        }
                        length_m = 0;
                        length_x += 1;
                        (in_m, in_x) = (false, true);
                    } else {
                        if in_x {
                            write!(out, "{length_x}X").expect("writing to a String cannot fail");
                        }
                        length_m += 1;
                        length_x = 0;
                        (in_m, in_x) = (true, false);
                    }
                    r += 1;
                    q += 1;
                }
            }
            b'I' => {
                q += i64::from(length);
                flush(&mut out, &mut in_m, &mut in_x, &mut length_m, &mut length_x);
                write!(out, "{length}I").expect("writing to a String cannot fail");
            }
            b'D' => {
                r += i64::from(length);
                flush(&mut out, &mut in_m, &mut in_x, &mut length_m, &mut length_x);
                write!(out, "{length}D").expect("writing to a String cannot fail");
            }
            _ => {}
        }
    }
    flush(&mut out, &mut in_m, &mut in_x, &mut length_m, &mut length_x);
    let end = query.len() as i32 - alignment.query_end - 1;
    if end > 0 {
        write!(out, "{end}S").expect("writing to a String cannot fail");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Opt-in: `SNV_SSW_CASES` holds `reference<TAB>query` lines and `SNV_SSW_EXPECTED`
    /// what the C library prints for them (`score ref_begin ref_end query_begin query_end
    /// cigar`, or `none`), e.g. from a harness linked against the original's
    /// `ssw.c`/`ssw_cpp.cpp`.
    #[test]
    fn matches_the_c_library() {
        let (Ok(cases), Ok(expected)) = (std::env::var("SNV_SSW_CASES"), std::env::var("SNV_SSW_EXPECTED")) else {
            return;
        };
        let cases = std::fs::read_to_string(cases).expect("cases");
        let expected = std::fs::read_to_string(expected).expect("expected");
        let mut checked = 0;
        for (case, want) in cases.lines().zip(expected.lines()) {
            let (reference, query) = case.split_once('\t').expect("two columns");
            let got = match Aligner::new(reference.as_bytes()).align(query.as_bytes()) {
                None => "none".to_owned(),
                Some(a) => format!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    a.score, a.ref_begin, a.ref_end, a.query_begin, a.query_end, a.cigar
                ),
            };
            assert_eq!(got, want, "case {case}");
            checked += 1;
        }
        assert!(checked > 0);
    }

    /// The SSE2 or NEON kernels (`sw_word`, `sw_byte` on x86-64 and aarch64) against the
    /// lane-by-lane emulation, both directions, with and without an early stop: on random
    /// references with `N`s, and on reads copied from them with substitutions and indels (long
    /// enough to overflow the 8-bit kernel) or unrelated.
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    #[test]
    fn vector_kernels_match_the_emulation() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let matrix = score_matrix();
        let bias = MISMATCH as u8;
        let base = |next: &mut dyn FnMut(usize) -> usize| if next(50) == 0 { 4 } else { next(4) as i8 };
        let mut checked = 0;
        for case in 0..4000 {
            let ref_len = 1 + next(400);
            let reference: Vec<i8> = (0..ref_len).map(|_| base(&mut next)).collect();
            let mut read: Vec<i8> = if case % 4 == 0 {
                (0..1 + next(300)).map(|_| base(&mut next)).collect()
            } else {
                let start = next(ref_len);
                reference[start..start + 1 + next((ref_len - start).min(300))].to_vec()
            };
            for _ in 0..next(10) {
                let at = next(read.len());
                match next(3) {
                    0 => read[at] = base(&mut next),
                    1 => read.insert(at, base(&mut next)),
                    _ if read.len() > 1 => {
                        read.remove(at);
                    }
                    _ => {}
                }
            }
            let read_len = read.len() as i32;
            let byte = profile(&read, &matrix, 16, |v| v.map_or(bias, |s| (i16::from(s) + i16::from(bias)) as u8));
            let word = profile(&read, &matrix, 8, |v| v.map_or(0, i16::from));
            let span = 1 + next(ref_len) as i32;
            for reverse in [false, true] {
                let full = lanes::sw_word(&reference, reverse, span, read_len, &word, u16::MAX);
                for terminate in [u16::MAX, full.score] {
                    let want = lanes::sw_word(&reference, reverse, span, read_len, &word, terminate);
                    let got = sw_word(&reference, reverse, span, read_len, &word, terminate);
                    assert_eq!((got.score, got.reference, got.read), (want.score, want.reference, want.read));
                    let terminate = terminate.min(255) as u8;
                    let want = lanes::sw_byte(&reference, reverse, span, read_len, &byte, terminate, bias);
                    let got = sw_byte(&reference, reverse, span, read_len, &byte, terminate, bias);
                    assert_eq!((got.score, got.reference, got.read), (want.score, want.reference, want.read));
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 16_000);
    }

    #[test]
    fn exact_match_is_all_equal() {
        let aligner = Aligner::new(b"TTTTACGTACGTGGCCAATTGG");
        let a = aligner.align(b"ACGTACGTGG").expect("aligned");
        assert_eq!((a.score, a.ref_begin, a.ref_end, a.query_begin, a.query_end), (40, 4, 13, 0, 9));
        assert_eq!(a.cigar, "10=");
    }

    #[test]
    fn mismatches_and_clips() {
        let aligner = Aligner::new(b"GGGGGACGTTCGTACGTAAAAA");
        let a = aligner.align(b"CCACGTACGTACGT").expect("aligned");
        assert_eq!(a.cigar, "2S4=1X7=");
        assert_eq!(a.score, 4 * 11 - 6);
    }
}
