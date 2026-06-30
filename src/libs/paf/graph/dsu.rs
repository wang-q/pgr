//! Union-find (union-by-rank + path compression) for graph node clustering.

/// Simple union-find structure used by [`super::PafGraph`] to compute the
/// transitive closure of aligned segments and to derive connected components
/// in the topology report.
pub(super) struct Dsu {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl Dsu {
    pub(super) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    pub(super) fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    pub(super) fn union(&mut self, a: usize, b: usize) {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        if self.rank[ra] == self.rank[rb] {
            self.rank[ra] += 1;
        }
    }
}
