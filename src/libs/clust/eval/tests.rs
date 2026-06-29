use super::*;
use crate::libs::pairmat::NamedMatrix;
use crate::libs::phylo::tree::Tree;

#[test]
fn test_tree_distance() {
    // ((A:2,B:4)g:2,(C:2,((D:3,E:1)h:1,F:2)i:1)j:1)k;
    // Construct tree manually or parse from string?
    // We don't have direct access to parser here easily unless we export it.
    // Tree::from_newick is likely available if io module supports it.
    // Tree::to_newick is available.
    // Let's assume we can parse it.
    // But Tree::from_file is available.
    // Let's create a dummy tree manually.

    let mut tree = Tree::new();
    let root = tree.add_node(); // k
    tree.set_root(root);

    let g = tree.add_node();
    let j = tree.add_node();
    tree.add_child(root, g).unwrap();
    tree.add_child(root, j).unwrap();

    let a = tree.add_node();
    let b = tree.add_node();
    tree.get_node_mut(a).unwrap().name = Some("A".to_string());
    tree.get_node_mut(b).unwrap().name = Some("B".to_string());
    tree.add_child(g, a).unwrap();
    tree.add_child(g, b).unwrap();

    // Set lengths
    // A:2, B:4. g:2 (to k?)
    // If we treat edges as parent->child length.
    tree.get_node_mut(a).unwrap().length = Some(2.0);
    tree.get_node_mut(b).unwrap().length = Some(4.0);
    tree.get_node_mut(g).unwrap().length = Some(2.0);

    let td = TreeDistance::new(tree);

    // Distance A-B = 2+4 = 6.
    assert_eq!(td.get_distance("A", "B"), 6.0);
}

#[test]
fn test_silhouette_score_simple() {
    // Data:
    // 0: 0.0 (C0)
    // 1: 1.0 (C1)
    // 2: 1.0 (C1)
    // 3: 2.0 (C1)
    // 4: 3.0 (C2)
    // 5: 3.0 (C2)

    let mut p = LabelMap::new();
    p.insert("0".to_string(), 0);
    p.insert("1".to_string(), 1);
    p.insert("2".to_string(), 1);
    p.insert("3".to_string(), 1);
    p.insert("4".to_string(), 2);
    p.insert("5".to_string(), 2);

    let names: Vec<String> = (0..6).map(|i| i.to_string()).collect();
    let mut dist_mat = NamedMatrix::new(names);
    let points: Vec<f32> = vec![0.0, 1.0, 1.0, 2.0, 3.0, 3.0];

    for i in 0..6 {
        for j in i + 1..6 {
            let d = (points[i] - points[j]).abs();
            let n1 = i.to_string();
            let n2 = j.to_string();
            dist_mat.set_by_name(&n1, &n2, d).unwrap();
        }
    }

    let score = silhouette_score(&p, &dist_mat);
    assert!((score - 0.5).abs() < 1e-6, "Score was {}", score);
}

#[test]
fn test_silhouette_score_single_cluster() {
    let mut p = LabelMap::new();
    p.insert("0".to_string(), 0);
    p.insert("1".to_string(), 0);

    let names = vec!["0".to_string(), "1".to_string()];
    let mut dist_mat = NamedMatrix::new(names);
    dist_mat.set_by_name("0", "1", 1.0).unwrap();

    let score = silhouette_score(&p, &dist_mat);
    assert_eq!(score, 0.0);
}

#[test]
fn test_silhouette_score_all_singletons() {
    // Sklearn behavior for all singletons is not strictly defined in docs but usually handled.
    // Our implementation returns 0.0 if n_clusters == n_samples
    let mut p = LabelMap::new();
    p.insert("0".to_string(), 0);
    p.insert("1".to_string(), 1);
    p.insert("2".to_string(), 2);

    let names = vec!["0".to_string(), "1".to_string(), "2".to_string()];
    let mut dist_mat = NamedMatrix::new(names);
    dist_mat.set_by_name("0", "1", 1.0).unwrap();
    dist_mat.set_by_name("0", "2", 1.0).unwrap();
    dist_mat.set_by_name("1", "2", 1.0).unwrap();

    let score = silhouette_score(&p, &dist_mat);
    assert_eq!(score, 0.0);
}

#[test]
fn test_davies_bouldin_score_simple() {
    // Cluster 1: A(0,0), B(0,1) -> Centroid (0, 0.5), Scatter = 0.5
    // Cluster 2: C(5,0), D(5,1) -> Centroid (5, 0.5), Scatter = 0.5
    // M12 = 5.0
    // R12 = (0.5+0.5)/5.0 = 0.2
    // DB = (0.2 + 0.2)/2 = 0.2

    let mut p = LabelMap::new();
    p.insert("A".to_string(), 1);
    p.insert("B".to_string(), 1);
    p.insert("C".to_string(), 2);
    p.insert("D".to_string(), 2);

    let mut data = HashMap::new();
    data.insert("A".to_string(), vec![0.0, 0.0]);
    data.insert("B".to_string(), vec![0.0, 1.0]);
    data.insert("C".to_string(), vec![5.0, 0.0]);
    data.insert("D".to_string(), vec![5.0, 1.0]);

    let coords = Coordinates { data, dim: 2 };

    let score = davies_bouldin_score(&p, &coords);
    assert!((score - 0.2).abs() < 1e-6, "Score was {}", score);
}

#[test]
fn test_evaluate_perfect() {
    let mut p1 = LabelMap::new();
    p1.insert("A".to_string(), 1);
    p1.insert("B".to_string(), 1);
    p1.insert("C".to_string(), 2);

    let mut p2 = LabelMap::new();
    p2.insert("A".to_string(), 10);
    p2.insert("B".to_string(), 10);
    p2.insert("C".to_string(), 20);

    let m = evaluate(&p1, &p2);
    assert_eq!(m.ari, 1.0);
    assert_eq!(m.ami, 1.0);
    assert_eq!(m.homogeneity, 1.0);
    assert_eq!(m.completeness, 1.0);
    assert_eq!(m.v_measure, 1.0);
    assert_eq!(m.fmi, 1.0);
    assert_eq!(m.nmi, 1.0);
}

#[test]
fn test_evaluate_disjoint() {
    // P1: {A,B}, {C,D} -> Labels: 1, 1, 2, 2
    // P2: {A,C}, {B,D} -> Labels: 1, 2, 1, 2
    // Contingency table is uniform:
    //      P2_1(AC) P2_2(BD)
    // P1_1(AB)  1(A)     1(B)
    // P1_2(CD)  1(C)     1(D)
    //
    // This is perfectly independent (orthogonal).
    // MI = 0.0
    // NMI = 0.0
    // ARI = -0.5 (Worse than random?) Let's check calculation:
    // sum_nij_2 = 0
    // sum_a_2 = 1 + 1 = 2
    // sum_b_2 = 1 + 1 = 2
    // n_2 = binom(4, 2) = 6
    // E[Index] = (2 * 2) / 6 = 4/6 = 0.666
    // Max[Index] = (2 + 2) / 2 = 2
    // Index = 0
    // ARI = (0 - 0.666) / (2 - 0.666) = -0.666 / 1.333 = -0.5
    // FMI = TP / sqrt(2 * 2) = 0 / 2 = 0.0

    let mut p1 = LabelMap::new();
    p1.insert("A".to_string(), 1);
    p1.insert("B".to_string(), 1);
    p1.insert("C".to_string(), 2);
    p1.insert("D".to_string(), 2);

    let mut p2 = LabelMap::new();
    p2.insert("A".to_string(), 1);
    p2.insert("C".to_string(), 1);
    p2.insert("B".to_string(), 2);
    p2.insert("D".to_string(), 2);

    let m = evaluate(&p1, &p2);
    assert!((m.ari + 0.5).abs() < 1e-6);
    assert!(m.mi.abs() < 1e-6);
    assert!(m.nmi.abs() < 1e-6);
    assert!(m.fmi.abs() < 1e-6);
}

#[test]
fn test_internal_indices_simple() {
    // Cluster 1: A(0,0), B(1,0) -> Centroid (0.5, 0)
    // Cluster 2: C(5,0), D(6,0) -> Centroid (5.5, 0)

    let mut p = LabelMap::new();
    p.insert("A".to_string(), 1);
    p.insert("B".to_string(), 1);
    p.insert("C".to_string(), 2);
    p.insert("D".to_string(), 2);

    let mut data = HashMap::new();
    data.insert("A".to_string(), vec![0.0, 0.0]);
    data.insert("B".to_string(), vec![1.0, 0.0]);
    data.insert("C".to_string(), vec![5.0, 0.0]);
    data.insert("D".to_string(), vec![6.0, 0.0]);

    let coords = Coordinates { data, dim: 2 };

    // Construct Distance Matrix for C-index
    let names = vec![
        "A".to_string(),
        "B".to_string(),
        "C".to_string(),
        "D".to_string(),
    ];
    let mut dist_mat = NamedMatrix::new(names);
    // Distances:
    // A-B: 1.0 (Intra)
    // C-D: 1.0 (Intra)
    // A-C: 5.0
    // A-D: 6.0
    // B-C: 4.0
    // B-D: 5.0
    // Sorted: 1, 1, 4, 5, 5, 6
    // N_W = 2 (A-B, C-D)
    // S_W = 1.0 + 1.0 = 2.0
    // S_min (sum of smallest 2) = 1.0 + 1.0 = 2.0
    // S_max (sum of largest 2) = 6.0 + 5.0 = 11.0
    // C-index = (2.0 - 2.0) / (11.0 - 2.0) = 0.0
    dist_mat.set_by_name("A", "B", 1.0).unwrap();
    dist_mat.set_by_name("C", "D", 1.0).unwrap();
    dist_mat.set_by_name("A", "C", 5.0).unwrap();
    dist_mat.set_by_name("A", "D", 6.0).unwrap();
    dist_mat.set_by_name("B", "C", 4.0).unwrap();
    dist_mat.set_by_name("B", "D", 5.0).unwrap();

    let c_index = c_index_score(&p, &dist_mat);
    assert_eq!(c_index, 0.0);

    // PBM:
    // Global Centroid: (12/4, 0) = (3, 0)
    // E_T: |0-3| + |1-3| + |5-3| + |6-3| = 3 + 2 + 2 + 3 = 10
    // E_W:
    //   C1: |0-0.5| + |1-0.5| = 0.5 + 0.5 = 1.0
    //   C2: |5-5.5| + |6-5.5| = 0.5 + 0.5 = 1.0
    //   Total E_W = 2.0
    // D_B: |0.5 - 5.5| = 5.0
    // K = 2
    // PBM = ( 1/2 * 10 / 2.0 * 5.0 )^2 = ( 0.5 * 5 * 5 )^2 = (12.5)^2 = 156.25
    let pbm = pbm_score(&p, &coords);
    assert!((pbm - 156.25).abs() < 1e-6, "PBM was {}", pbm);

    // Ball-Hall:
    // C1 mean dispersion: (|0-0.5|^2 + |1-0.5|^2) / 2 = (0.25+0.25)/2 = 0.25
    // C2 mean dispersion: (|5-5.5|^2 + |6-5.5|^2) / 2 = (0.25+0.25)/2 = 0.25
    // BH = (0.25 + 0.25) / 2 = 0.25
    let bh = ball_hall_score(&p, &coords);
    assert!((bh - 0.25).abs() < 1e-6, "Ball-Hall was {}", bh);

    // Xie-Beni:
    // WGSS = (0.25+0.25) + (0.25+0.25) = 1.0
    // min_sq_dist = (5.0)^2 = 25.0
    // N = 4
    // XB = 1.0 / (4 * 25.0) = 0.01
    let xb = xie_beni_score(&p, &coords);
    assert!((xb - 0.01).abs() < 1e-6, "Xie-Beni was {}", xb);

    // Wemmert-Gancarski:
    // C1 centroid G1=(0.5,0), C2 centroid G2=(5.5,0)
    // A(0,0): ||A-G1||=0.5, ||A-G2||=5.5. R(A)=0.5/5.5 = 1/11
    // B(1,0): ||B-G1||=0.5, ||B-G2||=4.5. R(B)=0.5/4.5 = 1/9
    // C(5,0): ||C-G2||=0.5, ||C-G1||=4.5. R(C)=0.5/4.5 = 1/9
    // D(6,0): ||D-G2||=0.5, ||D-G1||=5.5. R(D)=0.5/5.5 = 1/11
    // Mean R(C1) = (1/11 + 1/9)/2 = (9/99 + 11/99)/2 = 20/99/2 = 10/99
    // Mean R(C2) = 10/99
    // J1 = 1 - 10/99 = 89/99
    // J2 = 1 - 10/99 = 89/99
    // J = (2/4)*J1 + (2/4)*J2 = 0.5*J1 + 0.5*J2 = 89/99 = 0.898989...
    let wg = wemmert_gancarski_score(&p, &coords);
    let expected_wg = 89.0 / 99.0;
    assert!(
        (wg - expected_wg).abs() < 1e-6,
        "WG was {}, expected {}",
        wg,
        expected_wg
    );
}
