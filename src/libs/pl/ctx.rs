//! Pipeline working-directory context and path helpers.

use cmd_lib::*;
use intspan::absolute_path;

/// Resolve `path` to an absolute path string. `stdout` is passed through as-is.
pub fn abs_path_or_stdout(path: &str) -> anyhow::Result<String> {
    if path == "stdout" {
        Ok(path.to_string())
    } else {
        Ok(absolute_path(path)?.display().to_string())
    }
}

/// Shared pipeline context: current dir, pgr executable, and tempdir.
///
/// Created at the start of a pipeline; call [`PipelineCtx::enter`] to switch
/// into the tempdir and [`PipelineCtx::leave`] to restore the original cwd.
pub struct PipelineCtx {
    /// Original working directory, restored by `leave()`.
    pub curdir: std::path::PathBuf,
    /// Absolute path to the current `pgr` executable.
    pub pgr: String,
    /// Owned tempdir; dropped when the ctx is dropped.
    pub tempdir: tempfile::TempDir,
}

impl PipelineCtx {
    /// Create a new context with a tempdir using `prefix` (e.g. `"pgr_rm_"`).
    ///
    /// Prints the `==> Paths` info block.
    pub fn new(prefix: &str) -> anyhow::Result<Self> {
        let curdir = std::env::current_dir()?;
        let pgr = std::env::current_exe()?.display().to_string();
        let tempdir = tempfile::Builder::new().prefix(prefix).tempdir()?;
        let tempdir_str = tempdir.path().to_str().unwrap();

        run_cmd!(info "==> Paths")?;
        run_cmd!(info "    \"pgr\"     = ${pgr}")?;
        run_cmd!(info "    \"curdir\"  = ${curdir:?}")?;
        run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

        Ok(Self {
            curdir,
            pgr,
            tempdir,
        })
    }

    /// Resolve `p` to an absolute path string.
    pub fn abs_path(&self, p: &str) -> anyhow::Result<String> {
        Ok(absolute_path(p)?.display().to_string())
    }

    /// Switch the current working directory into the tempdir.
    pub fn enter(&self) -> anyhow::Result<()> {
        let tempdir_str = self.tempdir.path().to_str().unwrap();
        run_cmd!(info "==> Switch to tempdir")?;
        std::env::set_current_dir(tempdir_str)?;
        Ok(())
    }

    /// Restore the original working directory.
    pub fn leave(&self) -> anyhow::Result<()> {
        std::env::set_current_dir(&self.curdir)?;
        Ok(())
    }
}
