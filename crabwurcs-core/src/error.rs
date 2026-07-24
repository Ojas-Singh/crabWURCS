use thiserror::Error;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("failed to parse WURCS string at byte offset {offset}: {message}")]
    ParseError { offset: usize, message: String },

    #[error("WURCS string standardization failed: {0}")]
    StandardizationError(String),

    #[error("residue graph is malformed: {0} (e.g. dangling linkage, unreachable node from root)")]
    MalformedGraph(String),

    #[error("residue cannot be represented in standard WURCS: {0}")]
    UnrepresentableResidue(String),
}

pub type CoreResult<T> = Result<T, CoreError>;
