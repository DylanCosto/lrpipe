//! libstdc++'s `std::sort` (introsort, then insertion sort), which is not stable: equal
//! elements end up in the order this exact algorithm leaves them.

/// `_S_threshold`.
const THRESHOLD: usize = 16;

/// `std::sort(v.begin(), v.end(), less)`.
pub(crate) fn sort_by<T: Clone>(v: &mut [T], less: impl Fn(&T, &T) -> bool + Copy) {
    if v.len() < 2 {
        return;
    }
    let depth = 2 * v.len().ilog2() as usize;
    introsort_loop(v, depth, less);
    final_insertion_sort(v, less);
}

fn introsort_loop<T: Clone>(mut v: &mut [T], mut depth: usize, less: impl Fn(&T, &T) -> bool + Copy) {
    while v.len() > THRESHOLD {
        if depth == 0 {
            heap_sort(v, less);
            return;
        }
        depth -= 1;
        let cut = partition_pivot(v, less);
        let (left, right) = v.split_at_mut(cut);
        introsort_loop(right, depth, less);
        v = left;
    }
}

/// `__unguarded_partition_pivot`.
fn partition_pivot<T: Clone>(v: &mut [T], less: impl Fn(&T, &T) -> bool) -> usize {
    let mid = v.len() / 2;
    let (a, b, c) = (1, mid, v.len() - 1);
    // `__move_median_to_first(first, a, b, c)`.
    let pick = if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            b
        } else if less(&v[a], &v[c]) {
            c
        } else {
            a
        }
    } else if less(&v[a], &v[c]) {
        a
    } else if less(&v[b], &v[c]) {
        c
    } else {
        b
    };
    v.swap(0, pick);
    // `__unguarded_partition(first + 1, last, first)`.
    let (mut first, mut last) = (1, v.len());
    loop {
        while less(&v[first], &v[0]) {
            first += 1;
        }
        last -= 1;
        while less(&v[0], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn final_insertion_sort<T: Clone>(v: &mut [T], less: impl Fn(&T, &T) -> bool + Copy) {
    if v.len() > THRESHOLD {
        insertion_sort(&mut v[..THRESHOLD], less);
        for i in THRESHOLD..v.len() {
            linear_insert(v, i, less);
        }
    } else {
        insertion_sort(v, less);
    }
}

fn insertion_sort<T: Clone>(v: &mut [T], less: impl Fn(&T, &T) -> bool + Copy) {
    for i in 1..v.len() {
        if less(&v[i], &v[0]) {
            v[..=i].rotate_right(1);
        } else {
            linear_insert(v, i, less);
        }
    }
}

/// `__unguarded_linear_insert`.
fn linear_insert<T: Clone>(v: &mut [T], mut last: usize, less: impl Fn(&T, &T) -> bool) {
    let value = v[last].clone();
    while less(&value, &v[last - 1]) {
        v[last] = v[last - 1].clone();
        last -= 1;
    }
    v[last] = value;
}

/// `__partial_sort(first, last, last)`: `make_heap`, then `sort_heap`.
fn heap_sort<T: Clone>(v: &mut [T], less: impl Fn(&T, &T) -> bool + Copy) {
    let len = v.len();
    if len >= 2 {
        let mut parent = (len - 2) / 2;
        loop {
            let value = v[parent].clone();
            adjust_heap(v, parent, len, value, less);
            if parent == 0 {
                break;
            }
            parent -= 1;
        }
    }
    let mut last = len;
    while last > 1 {
        last -= 1;
        let value = v[last].clone();
        v[last] = v[0].clone();
        adjust_heap(&mut v[..last], 0, last, value, less);
    }
}

fn adjust_heap<T: Clone>(v: &mut [T], mut hole: usize, len: usize, value: T, less: impl Fn(&T, &T) -> bool) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if less(&v[second], &v[second - 1]) {
            second -= 1;
        }
        v[hole] = v[second].clone();
        hole = second;
    }
    if len & 1 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[hole] = v[second - 1].clone();
        hole = second - 1;
    }
    // `__push_heap`.
    let mut parent = hole.wrapping_sub(1) / 2;
    while hole > top && less(&v[parent], &value) {
        v[hole] = v[parent].clone();
        hole = parent;
        parent = hole.wrapping_sub(1) / 2;
    }
    v[hole] = value;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_by_key() {
        for n in [0usize, 1, 5, 16, 17, 40, 300] {
            let mut v: Vec<(i32, usize)> = (0..n).map(|i| (((i * 7919) % 13) as i32, i)).collect();
            sort_by(&mut v, |a, b| a.0 < b.0);
            assert!(v.windows(2).all(|w| w[0].0 <= w[1].0), "n = {n}");
        }
    }

    #[test]
    fn small_inputs_keep_ties_in_order() {
        let mut v = vec![(2, 0), (1, 1), (2, 2), (1, 3)];
        sort_by(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v, [(1, 1), (1, 3), (2, 0), (2, 2)]);
    }
}
