use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TreeError {
    /// Error during parsing (e.g., syntax error)
    ParseError {
        /// A human-readable message explaining the error
        message: String,
        /// The line number (1-based)
        line: usize,
        /// The column number (1-based)
        column: usize,
        /// The snippet of input where the error occurred
        snippet: String,
    },
    /// Logical error (e.g., cycle detected, invalid operation)
    LogicError(String),
}

impl fmt::Display for TreeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TreeError::ParseError {
                message,
                line,
                column,
                snippet,
            } => {
                write!(
                    f,
                    "Parse error at line {}, column {}:\n{}\nSnippet: \"{}\"",
                    line, column, message, snippet
                )
            }
            TreeError::LogicError(msg) => write!(f, "Tree logic error: {}", msg),
        }
    }
}

impl std::error::Error for TreeError {}
