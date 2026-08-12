//! External-crate API visibility check for the fq/asm migration (stage 1).
//!
//! After `fq`/`asm` move to anchr, the anchr crate depends on pgr for the
//! base layer (FASTA/FASTQ I/O, Phred encoding, k-mer, PAF, io/ds/loc/sys).
//! Integration tests compile against pgr as an external crate, so the `use`
//! statements below prove every symbol anchr needs is reachable through the
//! public API.

#![allow(dead_code, unused_imports)]

use pgr::libs::ds::radix_sort::radix_sort_bytes;
use pgr::libs::fmt::fq::{write_fa, write_fq};
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::pairs::PairReader;
use pgr::libs::fq::qual::{detect_quality_base, from_phred, to_phred, PHRED33, PHRED64};
use pgr::libs::kmer::count::count_keys;
use pgr::libs::kmer::key::Kmer;
use pgr::libs::kmer::qcheck::{check_read, CheckParams, ReadError};
use pgr::libs::kmer::quality::build_table;
use pgr::libs::kmer::{base_codes, canonical_keys, KmerTable};
use pgr::libs::nt::rev_comp;
use pgr::libs::par::ordered_map;
use pgr::libs::sys::mem_cap;

#[test]
fn migration_base_api_visible() {
    // Compile-time proof: every `use` above resolves, so the base layer is
    // available to an external crate (anchr).
    let _ = (PHRED33, PHRED64);
}
