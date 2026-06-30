use intspan::IntSpan;

/// Result of a Venn set-operation computation: exclusive element counts per set
/// and intersection counts ordered from lowest-order to highest-order intersections.
pub struct VennResult {
    /// Sizes of elements exclusive to each set (A only, B only, ...).
    pub excls: Vec<i32>,
    /// Sizes of intersections, ordered binary, then triple, ..., then the n-fold intersection.
    pub inter: Vec<i32>,
}

/// Compute Venn counts for 2 sets.
pub fn venn_sets_2(a: &IntSpan, b: &IntSpan) -> VennResult {
    let mut excls = Vec::new();
    let mut inter = Vec::new();

    // A ∩ B
    let i_ab = a.intersect(b).size();
    inter.push(i_ab);

    // A - B
    excls.push(a.diff(b).size());
    // B - A
    excls.push(b.diff(a).size());

    VennResult { excls, inter }
}

/// Compute Venn counts for 3 sets.
pub fn venn_sets_3(a: &IntSpan, b: &IntSpan, c: &IntSpan) -> VennResult {
    let mut excls = Vec::new();
    let mut inter = Vec::new();

    // A ∩ B ∩ C
    let i_abc = a.intersect(b).intersect(c);

    // Binary intersections minus triple intersection
    let sets_arr = [a, b, c];
    for i in 0..2 {
        for j in (i + 1)..=2 {
            let intersection = sets_arr[i].intersect(sets_arr[j]).diff(&i_abc).size();
            inter.push(intersection);
        }
    }

    inter.push(i_abc.size());

    // A - B - C
    excls.push(a.diff(b).diff(c).size());
    // B - A - C
    excls.push(b.diff(a).diff(c).size());
    // C - A - B
    excls.push(c.diff(a).diff(b).size());

    VennResult { excls, inter }
}

/// Compute Venn counts for 4 sets.
pub fn venn_sets_4(a: &IntSpan, b: &IntSpan, c: &IntSpan, d: &IntSpan) -> VennResult {
    let mut excls = Vec::new();
    let mut inter = Vec::new();

    // Quadruple intersection
    let i_abcd = a.intersect(b).intersect(c).intersect(d);

    // Binary intersections
    let sets_arr = [a, b, c, d];
    for i in 0..3 {
        for j in (i + 1)..=3 {
            let mut i_temp = sets_arr[i].intersect(sets_arr[j]);
            // Subtract all higher-order intersections containing these two sets
            for (k, _) in sets_arr.iter().enumerate() {
                if k != i && k != j {
                    i_temp.subtract(sets_arr[k]);
                }
            }
            inter.push(i_temp.size());
        }
    }

    // Triple intersections
    for i in 0..2 {
        for j in (i + 1)..3 {
            for k in (j + 1)..=3 {
                let i_temp = sets_arr[i]
                    .intersect(sets_arr[j])
                    .intersect(sets_arr[k])
                    .diff(&i_abcd);
                inter.push(i_temp.size());
            }
        }
    }

    // Quadruple intersection
    inter.push(i_abcd.size());

    // Exclusive elements
    excls.push(a.diff(b).diff(c).diff(d).size());
    excls.push(b.diff(a).diff(c).diff(d).size());
    excls.push(c.diff(a).diff(b).diff(d).size());
    excls.push(d.diff(a).diff(b).diff(c).size());

    VennResult { excls, inter }
}
