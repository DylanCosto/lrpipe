//! CPython's `list.sort` for short lists (fewer than 64 elements, where
//! timsort is a run detection plus binary insertion sort). Only `<` is used,
//! so values that are not ordered (NaN) end up where CPython puts them.

/// Sort `v` in place as CPython 3.12 does with the comparison `lt` (a < b).
pub fn sort_small<T>(v: &mut [T], lt: impl Fn(&T, &T) -> bool) {
    let n = v.len();
    assert!(n < 64, "only short lists are emulated");
    if n < 2 {
        return;
    }
    // count_run
    let mut run = 2;
    let descending = lt(&v[1], &v[0]);
    if descending {
        while run < n && lt(&v[run], &v[run - 1]) {
            run += 1;
        }
        v[..run].reverse();
    } else {
        while run < n && !lt(&v[run], &v[run - 1]) {
            run += 1;
        }
    }
    // binarysort of the rest
    for start in run..n {
        let (mut l, mut r) = (0, start);
        while l < r {
            let p = l + ((r - l) >> 1);
            if lt(&v[start], &v[p]) {
                r = p;
            } else {
                l = p + 1;
            }
        }
        v[l..=start].rotate_right(1);
    }
}

#[cfg(test)]
mod tests {
    use super::sort_small;

    #[test]
    fn like_cpython() {
        let mut v = [3.0, 1.0, 2.0];
        sort_small(&mut v, |a, b| a < b);
        assert_eq!(v, [1.0, 2.0, 3.0]);
        // sorted(..., key=lambda x: x[0]) in CPython 3.12
        let mut v = [(f64::NAN, 0), (1.0, 1), (0.5, 2)];
        sort_small(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v.map(|x| x.1), [0, 2, 1]);
        let mut v = [(1.0, 0), (f64::NAN, 1), (0.5, 2)];
        sort_small(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v.map(|x| x.1), [0, 1, 2]);
        // stable for equal keys
        let mut v = [(1, 'a'), (0, 'b'), (1, 'c'), (0, 'd')];
        sort_small(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v.map(|x| x.1), ['b', 'd', 'a', 'c']);
    }
}
