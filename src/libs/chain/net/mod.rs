//! Chain-to-Net construction logic and Net format I/O.
//!
//! Builds a Net (hierarchical gap-fill tree) from a set of Chains by inserting
//! each chain's alignment blocks into the target/query chromosome gap trees.
//! Also provides readers/writers for the UCSC Net text format.
//!
//! Module layout:
//! * [`types`] — Net data structures (NetNode, Space, Gap, Fill, Chrom).
//! * [`builder`] — ChainNet builder (inserts chains into gap trees).
//! * [`reader`] — UCSC Net text format reader.
//! * [`writer`] — UCSC Net text format writer (filtered, with subchain scoring).
//! * [`finalize`] — sort + recompute o_start/o_end from chain data.
//! * [`syntenic`] — `classify_syntenic` for query-side duplication depth classification.

pub mod builder;
pub mod class;
pub mod filter;
pub mod finalize;
pub mod reader;
pub mod subset;
pub mod syntenic;
pub mod to_axt;
pub mod types;
pub mod writer;

pub use builder::ChainNet;
pub use class::{collect_stats_fill, collect_stats_gap, Stats};
pub use filter::{filter_chrom, prune_gap, FilterCriteria};
pub use finalize::finalize_net;
pub use reader::read_nets;
pub use subset::{subset_nets, SubsetOptions};
pub use syntenic::classify_syntenic;
pub use to_axt::net_to_axt;
pub use types::{Chrom, Fill, Gap, NetNode, Space};
pub use writer::{range_intersection, write_net, write_net_file, write_sorted_net};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ucsc_net_format_compatibility() {
        let net_data = "\
net chr2L 23011544
 fill 6004 3278 chrXR_group3a - 1396397 2164 id 25606 score 23114 ali 782 qDup 576 type top tN 0 qN 0 tR 36 qR 0 tTrf 0 qTrf 0
  gap 6065 2 chrXR_group3a - 1398498 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 6096 1485 chrXR_group3a - 1397572 897 tN 0 qN 0 tR 36 qR 0 tTrf 0 qTrf 0
   fill 6096 513 chrU - 5570675 533 id 48675 score 4435 ali 465 qDup 533 type nonSyn tN 0 qN 0 tR 0 qR 13 tTrf 0 qTrf 0
    gap 6116 8 chrU - 5571188 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 6156 5 chrU - 5571156 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 6184 3 chrU - 5571133 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 6212 18 chrU - 5571106 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 6244 9 chrU - 5571092 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 6340 2 chrU - 5570996 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 6515 3 chrU - 5570771 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 7623 1 chrXR_group3a - 1397530 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 7664 1007 chrXR_group3a - 1397008 482 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
   fill 7664 382 chrXL_group1e - 8262003 506 id 25608 score 10609 ali 364 qDup 506 type nonSyn tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 7784 4 chrXL_group1e - 8262361 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 7792 3 chrXL_group1e - 8262357 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 7921 2 chrXL_group1e - 8262126 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
    gap 7949 9 chrXL_group1e - 8262092 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 8693 1 chrXR_group3a - 1396985 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
 fill 9833 1251 chrU - 5562980 1239 id 48675 score 10720 ali 1124 qDup 1094 type top tN 0 qN 0 tR 22 qR 88 tTrf 0 qTrf 0
  gap 9966 7 chrU - 5564075 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 10015 3 chrU - 5564030 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 10088 2 chrU - 5563957 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
  gap 10101 8 chrU - 5563946 0 tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
";
        let reader = std::io::Cursor::new(net_data);
        let chroms = read_nets(reader).unwrap();

        assert_eq!(chroms.len(), 1);
        let chrom = &chroms[0];
        assert_eq!(chrom.name, "chr2L");
        assert_eq!(chrom.size, 23011544);

        let mut writer = Vec::new();
        chrom.write(&mut writer).unwrap();
        let output = String::from_utf8(writer).unwrap();

        assert_eq!(output, net_data);
    }
}
