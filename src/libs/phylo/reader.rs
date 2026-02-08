use super::tree::Tree;
use std::io::Read;

/// Read a Newick tree from a file.
///
/// # Arguments
/// * `infile` - Path to the input file (or "-" for stdin).
///
/// # Example
/// ```
/// // usage in CLI:
/// // let tree = pgr::libs::phylo::reader::from_file("path/to/tree.nwk");
/// ```
pub fn from_file(infile: &str) -> Tree {
    let mut reader = intspan::reader(infile);
    let mut newick = String::new();
    reader.read_to_string(&mut newick).expect("Read error");
    Tree::from_newick(newick.as_str()).unwrap()
}
