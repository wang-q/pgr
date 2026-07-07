//! Pipeline working-directory context and path helpers.

use cmd_lib::run_cmd;
use intspan::absolute_path;

/// Resolve `path` to an absolute path string. `stdout` is passed through as-is.
pub fn abs_path_or_stdout(path: &str) -> anyhow::Result<String> {
    if path == "stdout" {
        Ok(path.to_string())
    } else {
        Ok(absolute_path(path)?.display().to_string())
    }
}

/// Shared pipeline context: pgr executable and tempdir.
///
/// Created at the start of a pipeline; call [`PipelineCtx::enter`] to switch
/// into the tempdir. The returned [`super::CwdGuard`] restores the original
/// working directory on drop, so CWD is always restored — even on error.
pub struct PipelineCtx {
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
        let pgr = crate::libs::io::current_exe_string()?;
        let tempdir = tempfile::Builder::new().prefix(prefix).tempdir()?;
        let tempdir_str = tempdir.path().to_str().unwrap();

        run_cmd!(info "==> Paths")?;
        run_cmd!(info "    \"pgr\"     = ${pgr}")?;
        run_cmd!(info "    \"curdir\"  = ${curdir:?}")?;
        run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

        Ok(Self { pgr, tempdir })
    }

    /// Resolve `p` to an absolute path string.
    pub fn abs_path(&self, p: &str) -> anyhow::Result<String> {
        Ok(absolute_path(p)?.display().to_string())
    }

    /// Switch the current working directory into the tempdir.
    ///
    /// Returns a [`super::CwdGuard`] whose `Drop` restores the previous
    /// working directory, ensuring cleanup even when the pipeline errors out.
    pub fn enter(&self) -> anyhow::Result<super::CwdGuard> {
        let tempdir_str = self.tempdir.path().to_str().unwrap();
        run_cmd!(info "==> Switch to tempdir")?;
        super::CwdGuard::enter(tempdir_str)
    }
}
