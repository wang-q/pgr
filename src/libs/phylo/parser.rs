use super::error::TreeError;
use super::node::NodeId;
use super::tree::Tree;
use nom::{
    branch::alt,
    bytes::complete::{is_not, take_while},
    character::complete::{char, digit1, multispace0},
    combinator::{cut, map, map_res, opt, recognize},
    error::{context, ContextError, ErrorKind, FromExternalError, ParseError},
    multi::{many1, separated_list1},
    sequence::{delimited, preceded},
    IResult, Offset, Parser,
};
use std::collections::BTreeMap;

// ================================================================================================
// Error Handling Structures
// ================================================================================================

#[derive(Clone, Debug, PartialEq)]
pub enum DetailedErrorKind {
    Context(&'static str),
    Nom(ErrorKind),
}

/// A custom error type for nom that accumulates context and error kinds.
/// This allows for more informative error messages when parsing fails.
#[derive(Clone, Debug, PartialEq)]
pub struct DetailedError<'a> {
    pub errors: Vec<(&'a str, DetailedErrorKind)>,
}

impl<'a> ParseError<&'a str> for DetailedError<'a> {
    fn from_error_kind(input: &'a str, kind: ErrorKind) -> Self {
        DetailedError {
            errors: vec![(input, DetailedErrorKind::Nom(kind))],
        }
    }

    fn append(input: &'a str, kind: ErrorKind, mut other: Self) -> Self {
        other.errors.push((input, DetailedErrorKind::Nom(kind)));
        other
    }
}

impl<'a> ContextError<&'a str> for DetailedError<'a> {
    fn add_context(input: &'a str, ctx: &'static str, mut other: Self) -> Self {
        other.errors.push((input, DetailedErrorKind::Context(ctx)));
        other
    }
}

impl<'a, E> FromExternalError<&'a str, E> for DetailedError<'a> {
    fn from_external_error(input: &'a str, kind: ErrorKind, _e: E) -> Self {
        DetailedError {
            errors: vec![(input, DetailedErrorKind::Nom(kind))],
        }
    }
}

// ================================================================================================
// Intermediate Structure
// ================================================================================================

/// `ParsedNode` is a temporary recursive structure used during parsing.
/// It mirrors the structure of a Newick tree node but exists independently of the final `Tree` arena.
///
/// Parsing a recursive structure like Newick is easier when building a recursive data type.
/// However, the final `Tree` structure in this library uses an arena-based (flat vector) approach
/// for better performance and memory locality.
///
/// After parsing is complete, this structure is converted into the flat, arena-based `Tree`
/// via the `to_tree` method.
#[derive(Debug)]
struct ParsedNode {
    name: Option<String>,
    length: Option<f64>,
    properties: Option<BTreeMap<String, String>>, // For NHX comments: [&&NHX:key=value]
    children: Vec<ParsedNode>,
}

impl ParsedNode {
    fn new() -> Self {
        Self {
            name: None,
            length: None,
            properties: None,
            children: Vec::new(),
        }
    }

    /// Converts this recursive `ParsedNode` into nodes in the provided `Tree` arena.
    /// Returns the `NodeId` of the created node in the arena.
    ///
    /// This function recursively traverses the `ParsedNode` tree, creating `Node`s in the `Tree` struct
    /// and linking them together.
    fn to_tree(self, tree: &mut Tree) -> NodeId {
        let id = tree.add_node();
        for child in self.children {
            let child_id = child.to_tree(tree);
            // The unwrap here is safe because `id` was just created and exists in the tree.
            tree.add_child(id, child_id).unwrap();
        }
        if let Some(node) = tree.get_node_mut(id) {
            node.name = self.name;
            node.length = self.length;
            node.properties = self.properties;
        }
        id
    }
}

// ================================================================================================
// Parsers
// ================================================================================================

// 1. Whitespace eater
// In nom 8, we return `impl Parser` instead of `impl FnMut`.
// This parser wraps another parser and ignores surrounding whitespace (spaces, tabs, newlines).
// It's used extensively to make the parser robust against formatting variations.
fn ws<'a, F, O, E>(inner: F) -> impl Parser<&'a str, Output = O, Error = E>
where
    F: Parser<&'a str, Output = O, Error = E>,
    E: ParseError<&'a str>,
{
    delimited(multispace0, inner, multispace0)
}

// 2. Label
// Parses a node label/name.
// Supports:
// - Unquoted strings (stops at reserved chars: "():;,[]")
// - Single quoted strings ('example name') - internal single quotes can be escaped as ''
// - Double quoted strings ("example name") - internal double quotes can be escaped as ""
fn parse_label(input: &str) -> IResult<&str, String, DetailedError<'_>> {
    // Unquoted labels cannot contain Newick structural characters
    let unquoted = map(
        // Take characters until a reserved Newick character is found
        take_while(|c: char| !"():;,[]".contains(c)),
        |s: &str| s.trim().to_string(),
    );

    // Single quoted labels: 'Homo sapiens'
    // Two single quotes inside represent one single quote: 'O''Brien' -> O'Brien
    let single_quoted = delimited(
        char('\''),
        map(is_not("'"), |s: &str| s.replace("''", "'")),
        char('\''),
    );

    // Double quoted labels: "Homo sapiens"
    // Two double quotes inside represent one double quote: "He said ""Hello""" -> He said "Hello"
    let double_quoted = delimited(
        char('"'),
        map(is_not("\""), |s: &str| s.replace("\"\"", "\"")),
        char('"'),
    );

    // Try quoted formats first, then fall back to unquoted
    context("label", alt((single_quoted, double_quoted, unquoted))).parse(input)
}

// 3. Length
// Parses the branch length, which follows a colon (e.g., ":0.123").
// Supports standard floating point formats including scientific notation.
fn parse_length(input: &str) -> IResult<&str, f64, DetailedError<'_>> {
    context(
        "length",
        preceded(
            ws(char(':')), // Lengths must start with ':'
            // Use `cut` to prevent backtracking if we found a ':' but failed to parse the number.
            // This gives a better error message ("expected float" instead of trying other branches).
            cut(map_res(
                recognize((
                    opt(char('-')),
                    digit1,
                    opt((char('.'), digit1)),
                    opt((
                        alt((char('e'), char('E'))),
                        opt(alt((char('+'), char('-')))),
                        digit1,
                    )),
                )),
                |s: &str| s.parse::<f64>(),
            )),
        ),
    )
    .parse(input)
}

// 4. Comment
// Parses Newick comments enclosed in square brackets: [comment].
// Specifically handles NHX (New Hampshire eXtended) format comments: [&&NHX:key=value:...]
// Returns:
// - Some(map) if it's an NHX comment with properties
// - None if it's a regular comment (ignored)
fn parse_comment(
    input: &str,
) -> IResult<&str, Option<BTreeMap<String, String>>, DetailedError<'_>> {
    let comment_content = delimited(ws(char('[')), is_not("]"), char(']'));

    context(
        "comment",
        map(opt(comment_content), |content: Option<&str>| {
            if let Some(s) = content {
                // Check for NHX format signature
                // Example: [&&NHX:S=human:E=1.5]
                if s.starts_with("&&NHX") {
                    let mut props = BTreeMap::new();
                    for part in s.split(':') {
                        if part == "&&NHX" {
                            continue;
                        }
                        if let Some((k, v)) = part.split_once('=') {
                            props.insert(k.to_string(), v.to_string());
                        }
                    }
                    if !props.is_empty() {
                        return Some(props);
                    }
                } else {
                    // Try to parse simple Key=Value properties
                    // Example: [S=Gorilla]
                    let mut props = BTreeMap::new();
                    for part in s.split_whitespace() {
                        if let Some((k, v)) = part.split_once('=') {
                            props.insert(k.to_string(), v.to_string());
                        }
                    }
                    if !props.is_empty() {
                        return Some(props);
                    }
                }
            }
            None
        }),
    )
    .parse(input)
}

// 5. Subtree
// Recursive parser for a tree node and its children.
// General Newick Structure: (child1, child2, ...)Label:Length[Comment]
fn parse_subtree(input: &str) -> IResult<&str, ParsedNode, DetailedError<'_>> {
    // 1. Children: optional list of subtrees enclosed in parens
    // Example: (A:0.1, B:0.2)
    let (input, children) = context(
        "children",
        opt(delimited(
            ws(char('(')),
            separated_list1(ws(char(',')), parse_subtree), // Comma-separated list of subtrees
            ws(char(')')),
        )),
    )
    .parse(input)?;

    // 2. Label: optional node name
    // Example: Homo_sapiens
    let (input, label) = opt(parse_label).parse(input)?;

    // 3. Properties/Comments/Length:
    // Newick allows comments before or after length, so we parse both.
    // Example: :0.1[&&NHX:...] or [&&NHX:..]:0.1
    let (input, comment1) = parse_comment(input)?;
    let (input, length) = opt(parse_length).parse(input)?;
    let (input, comment2) = parse_comment(input)?;

    // Construct the intermediate ParsedNode
    let mut node = ParsedNode::new();
    if let Some(c) = children {
        node.children = c;
    }
    if let Some(l) = label {
        if !l.is_empty() {
            node.name = Some(l);
        }
    }
    node.length = length;

    // Merge properties from comments found before and after length
    if comment1.is_some() || comment2.is_some() {
        let mut props = BTreeMap::new();
        if let Some(p) = comment1 {
            props.extend(p);
        }
        if let Some(p) = comment2 {
            props.extend(p);
        }
        node.properties = Some(props);
    }

    Ok((input, node))
}

// ================================================================================================
// Entry Points
// ================================================================================================

// 6. Entry point - Single Tree
/// Parses a single Newick tree string.
/// Expects the tree to end with a semicolon ';'.
///
/// # Arguments
/// * `input` - The Newick string to parse.
///
/// # Returns
/// * `Result<Tree, TreeError>` - The parsed tree or an error.
pub fn parse_newick(input: &str) -> Result<Tree, TreeError> {
    let mut parser = (ws(parse_subtree), ws(char(';')));

    match parser.parse(input) {
        Ok((_, (root_node, _))) => {
            let mut tree = Tree::new();
            let root_id = root_node.to_tree(&mut tree);
            tree.set_root(root_id);
            Ok(tree)
        }
        Err(nom::Err::Error(e)) | Err(nom::Err::Failure(e)) => Err(make_tree_error(input, e)),
        Err(nom::Err::Incomplete(_)) => Err(TreeError::ParseError {
            message: "Incomplete input".to_string(),
            line: 0,
            column: 0,
            snippet: "".to_string(),
        }),
    }
}

// Entry point - Multiple Trees
/// Parses a string containing multiple Newick trees.
/// Handles standard trees ending in ';' as well as "garbage" or comments
/// (like file headers) that are enclosed in square brackets `[...]` but are not part of a tree.
///
/// This is useful for parsing files that contain multiple trees or have metadata headers.
pub fn parse_newick_multi(input: &str) -> Result<Vec<Tree>, TreeError> {
    // A valid tree is a subtree followed by a semicolon
    let valid_tree = map((ws(parse_subtree), ws(char(';'))), |(root, _)| Some(root));

    // "Garbage" blocks are top-level comments [ ... ] that are ignored.
    // Some tree files (like from Nexus) might have headers in comments.
    let garbage = map(
        ws(delimited(char('['), take_while(|c| c != ']'), char(']'))),
        |_| None,
    );

    // Parse many occurrences of either valid trees or garbage
    let mut parser = many1(alt((valid_tree, garbage)));

    match parser.parse(input) {
        Ok((_, trees_data)) => {
            let mut trees = Vec::new();
            for root_opt in trees_data {
                if let Some(root_node) = root_opt {
                    let mut tree = Tree::new();
                    let root_id = root_node.to_tree(&mut tree);
                    tree.set_root(root_id);
                    trees.push(tree);
                }
            }
            Ok(trees)
        }
        Err(nom::Err::Error(e)) | Err(nom::Err::Failure(e)) => Err(make_tree_error(input, e)),
        Err(nom::Err::Incomplete(_)) => Err(TreeError::ParseError {
            message: "Incomplete input".to_string(),
            line: 0,
            column: 0,
            snippet: "".to_string(),
        }),
    }
}

// Helper to convert nom errors into friendly TreeError
fn make_tree_error(input: &str, e: DetailedError) -> TreeError {
    let (remaining, _) = e.errors.first().unwrap();
    let offset = input.offset(remaining);

    // Calculate line/col
    let prefix = &input[..offset];
    let line = prefix.chars().filter(|&c| c == '\n').count() + 1;
    let last_newline = prefix.rfind('\n').map(|p| p + 1).unwrap_or(0);
    let column = offset - last_newline + 1;

    let mut msg = String::new();
    for (_, kind) in e.errors.iter().rev() {
        match kind {
            DetailedErrorKind::Context(ctx) => {
                msg.push_str(&format!("while parsing {}:\n", ctx));
            }
            DetailedErrorKind::Nom(k) => {
                msg.push_str(&format!("  error: {:?}\n", k));
            }
        }
    }

    TreeError::ParseError {
        message: msg,
        line,
        column,
        snippet: remaining.chars().take(50).collect(),
    }
}

impl Tree {
    /// Parse a Newick string into a Tree.
    ///
    /// # Example
    /// ```
    /// use pgr::libs::phylo::tree::Tree;
    ///
    /// // Successful parse
    /// let input = "(A:0.1,B:0.2)Root;";
    /// let tree = Tree::from_newick(input).unwrap();
    /// assert_eq!(tree.len(), 3);
    ///
    /// // Error handling
    /// let invalid_input = "(A,B:invalid)C;";
    /// let result = Tree::from_newick(invalid_input);
    /// assert!(result.is_err());
    /// println!("Error: {}", result.err().unwrap());
    /// ```
    pub fn from_newick(input: &str) -> Result<Self, TreeError> {
        parse_newick(input)
    }

    pub fn from_newick_multi(input: &str) -> Result<Vec<Self>, TreeError> {
        parse_newick_multi(input)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parser_simple() {
        let input = "(A,B)C;";
        let tree = Tree::from_newick(input).unwrap();
        assert_eq!(tree.len(), 3);

        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        assert_eq!(root.name.as_deref(), Some("C"));
        assert_eq!(root.children.len(), 2);
    }

    #[test]
    fn test_parser_lengths() {
        let input = "(A:0.1, B:0.2e-1)Root:100;";
        let tree = Tree::from_newick(input).unwrap();

        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        assert_eq!(root.name.as_deref(), Some("Root"));
        assert_eq!(root.length, Some(100.0));

        let child1 = tree.get_node(root.children[0]).unwrap();
        assert_eq!(child1.name.as_deref(), Some("A"));
        assert_eq!(child1.length, Some(0.1));

        let child2 = tree.get_node(root.children[1]).unwrap();
        assert_eq!(child2.name.as_deref(), Some("B"));
        assert_eq!(child2.length, Some(0.02)); // 0.2e-1
    }

    #[test]
    fn test_parser_nhx() {
        let input = "(A:0.1,B:0.2)n1[&&NHX:S=human:E=1.5];";
        let tree = Tree::from_newick(input).unwrap();

        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        assert_eq!(root.name.as_deref(), Some("n1"));

        let props = root.properties.as_ref().unwrap();
        assert_eq!(props.get("S").map(|s| s.as_str()), Some("human"));
        assert_eq!(props.get("E").map(|s| s.as_str()), Some("1.5"));
    }

    #[test]
    fn test_parser_whitespace() {
        let input = "  (  A : 0.1 ,  B  )  ;  ";
        let tree = Tree::from_newick(input).unwrap();
        assert_eq!(tree.len(), 3);
    }

    #[test]
    fn test_parser_multiline_whitespace() {
        // User request: support extensive whitespace and newlines for readability
        let input = "
        (
            A : 0.1,
            B : 0.2
        ) Root ;
        ";
        let tree = Tree::from_newick(input).unwrap();
        assert_eq!(tree.len(), 3);

        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        assert_eq!(root.name.as_deref(), Some("Root"));

        let c0 = tree.get_node(root.children[0]).unwrap();
        assert_eq!(c0.name.as_deref(), Some("A"));
        assert_eq!(c0.length, Some(0.1));

        let c1 = tree.get_node(root.children[1]).unwrap();
        assert_eq!(c1.name.as_deref(), Some("B"));
        assert_eq!(c1.length, Some(0.2));
    }

    #[test]
    fn test_parser_complex_formatting() {
        // A more complex example with nested structure and comments across lines
        // Note: Commas must come AFTER the node info (label:length[comment]),
        // so [Comment] must be before the comma if it belongs to that node.
        let input = "
        (
            (
                'Human' : 0.1    [Comment on Human],
                'Chimp' : 0.12    [Comment on Chimp]
            )Hominidae : 0.5,
            'Gorilla' : 0.6
        )Hominoidea;
        ";
        let tree = Tree::from_newick(input).unwrap();
        // Structure:
        //        Hominoidea
        //       /          \
        //  Hominidae      Gorilla
        //    /   \
        // Human Chimp

        assert_eq!(tree.len(), 5);
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        assert_eq!(root.name.as_deref(), Some("Hominoidea"));

        let gorilla = tree.get_node(root.children[1]).unwrap();
        assert_eq!(gorilla.name.as_deref(), Some("Gorilla"));

        let hominidae = tree.get_node(root.children[0]).unwrap();
        assert_eq!(hominidae.name.as_deref(), Some("Hominidae"));
        assert_eq!(hominidae.children.len(), 2);
    }

    #[test]
    fn test_parser_whitespace_details() {
        // Case 1: Spaces around unquoted labels
        // "( A , B )" -> names should be "A" and "B" without spaces
        let input = "( A , B )Root;";
        let tree = Tree::from_newick(input).unwrap();
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        let c0 = tree.get_node(root.children[0]).unwrap();
        assert_eq!(c0.name.as_deref(), Some("A"));
        let c1 = tree.get_node(root.children[1]).unwrap();
        assert_eq!(c1.name.as_deref(), Some("B"));

        // Case 2: Spaces inside quoted labels (should be preserved by parser)
        let input = "(' A ', ' B ')Root;";
        let tree = Tree::from_newick(input).unwrap();
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        let c0 = tree.get_node(root.children[0]).unwrap();
        let c1 = tree.get_node(root.children[1]).unwrap();

        assert_eq!(c0.name.as_deref(), Some(" A "));
        assert_eq!(c1.name.as_deref(), Some(" B "));

        // Case 3: Mixed spaces
        let input = "(  A:0.1  ,  B:0.2  )  Root  ;";
        let tree = Tree::from_newick(input).unwrap();
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        let c0 = tree.get_node(root.children[0]).unwrap();
        let c1 = tree.get_node(root.children[1]).unwrap();

        assert_eq!(c0.name.as_deref(), Some("A"));
        assert_eq!(c1.name.as_deref(), Some("B"));
        assert_eq!(root.name.as_deref(), Some("Root"));
    }

    #[test]
    fn test_parser_user_scenario_emulation() {
        // Emulate the user's cleaning logic to show it's mostly redundant for unquoted,
        // but works for quoted if that's what they want.
        let input = "(' A ', B )Root;";
        let mut tree = Tree::from_newick(input).unwrap();

        // Before cleaning
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        let c0 = tree.get_node(root.children[0]).unwrap(); // ' A '
        let c1 = tree.get_node(root.children[1]).unwrap(); // B

        assert_eq!(c0.name.as_deref(), Some(" A "));
        assert_eq!(c1.name.as_deref(), Some("B"));

        // User's cleaning logic
        let root_id = tree.get_root().unwrap();
        let traversal = tree.preorder(&root_id).unwrap();

        for id in traversal {
            let node = tree.get_node_mut(id).unwrap();
            if let Some(ref name) = node.name.clone() {
                let trimmed = name.trim().to_string();
                if trimmed.is_empty() {
                    node.name = None;
                } else {
                    node.set_name(trimmed);
                }
            }
        }

        // After cleaning
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();
        let c0 = tree.get_node(root.children[0]).unwrap();
        let c1 = tree.get_node(root.children[1]).unwrap();

        assert_eq!(c0.name.as_deref(), Some("A")); // Trimmed!
        assert_eq!(c1.name.as_deref(), Some("B")); // Unchanged
    }

    #[test]
    fn test_parser_double_quoted() {
        let input = "(\"Homo sapiens\":0.1, \"Mus musculus\":0.2);";
        let tree = Tree::from_newick(input).unwrap();
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();

        let c1 = tree.get_node(root.children[0]).unwrap();
        assert_eq!(c1.name.as_deref(), Some("Homo sapiens"));

        let c2 = tree.get_node(root.children[1]).unwrap();
        assert_eq!(c2.name.as_deref(), Some("Mus musculus"));
    }

    #[test]
    fn test_parser_quoted() {
        let input = "('Homo sapiens':0.1, 'Mus musculus':0.2);";
        let tree = Tree::from_newick(input).unwrap();
        let root = tree.get_node(tree.get_root().unwrap()).unwrap();

        let c1 = tree.get_node(root.children[0]).unwrap();
        assert_eq!(c1.name.as_deref(), Some("Homo sapiens"));

        let c2 = tree.get_node(root.children[1]).unwrap();
        assert_eq!(c2.name.as_deref(), Some("Mus musculus"));
    }

    #[test]
    fn test_parser_error() {
        // Case 1: Missing semicolon
        let input = "(A,B)C";
        let res = Tree::from_newick(input);
        match res {
            Err(TreeError::ParseError { line, column, .. }) => {
                assert_eq!(line, 1);
                // (A,B)C -> length 6. Expects ; at col 7.
                assert_eq!(column, 7);
            }
            _ => panic!("Expected ParseError, got {:?}", res),
        }

        // Case 2: Invalid length
        let input2 = "(A,B:invalid)C;";
        let res2 = Tree::from_newick(input2);
        match res2 {
            Err(TreeError::ParseError { line, message, .. }) => {
                assert_eq!(line, 1);
                assert!(message.contains("length"));
            }
            _ => panic!("Expected ParseError, got {:?}", res2),
        }
    }
}
