//! Stage-cost profile for super-mer counting (task 1.1, stage 0).
//!
//! Prints per-k wall time and, with `PGR_SUPERMER_TIMING=1`, the
//! stage-1 pack / stage-1 sort / stage-2 expand / stage-2 sort breakdown
//! emitted by `kmer::supermer`. Defaults to the MG1655 20x simulated read
//! set (150 bp windows, 7 bp stride over 1 Mb); pass a FASTA path as the
//! first argument to profile real reads instead.

use pgr::libs::kmer::supermer;
use pgr::libs::pgi::build::read_fasta;

fn main() -> anyhow::Result<()> {
    let (seqs, label) = match std::env::args().nth(1) {
        Some(path) => {
            let seqs: Vec<Vec<u8>> = read_fasta(&path)?.into_iter().map(|(_, s)| s).collect();
            (seqs, path)
        }
        None => {
            let genome = read_fasta(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/genome/mg1655.fa.gz"
            ))?;
            let contig = genome
                .iter()
                .map(|(_, s)| s.as_slice())
                .find(|s| s.len() >= 1_000_000)
                .expect("mg1655 has a >= 1 Mb contig");
            let reads: Vec<Vec<u8>> = (0..1_000_000usize - 150)
                .step_by(7)
                .map(|p| contig[p..p + 150].to_vec())
                .collect();
            (reads, "mg1655 reads20x".to_string())
        }
    };
    let bases: usize = seqs.iter().map(Vec::len).sum();
    eprintln!("data: {label}, seqs={}, bases={bases}", seqs.len());
    for k in [21usize, 41, 61, 81] {
        let t0 = std::time::Instant::now();
        let table = supermer::build_table(&seqs, k)?;
        let wall = t0.elapsed();
        eprintln!(
            "k={k}: wall={:.3}s unique={}",
            wall.as_secs_f64(),
            table.counts.len()
        );
    }
    Ok(())
}
