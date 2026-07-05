//! Shared helpers for `pgr pl` pipeline subcommands.
//!
//! Pure pipeline orchestration logic (no clap dependency): workflow context,
//! path resolution, and external-tool driver functions (FastK / Profex / spanr).

mod ctx;
mod repeat;

pub use ctx::{abs_path_or_stdout, PipelineCtx};
pub use repeat::{
    parse_trf_output, run_profex_per_chr, run_repeat_pipeline, run_repeat_spanr_pipeline,
    RepeatOpts,
};

use std::path::PathBuf;

/// RAII guard that restores the working directory on drop.
pub struct CwdGuard {
    prev_dir: PathBuf,
}

impl CwdGuard {
    /// Change to `new_dir` and return a guard that restores the previous
    /// directory on drop.
    pub fn enter(new_dir: &str) -> anyhow::Result<Self> {
        let prev_dir = std::env::current_dir()?;
        std::env::set_current_dir(new_dir)?;
        Ok(Self { prev_dir })
    }
}

impl Drop for CwdGuard {
    fn drop(&mut self) {
        if let Err(e) = std::env::set_current_dir(&self.prev_dir) {
            log::warn!("failed to restore working directory: {}", e);
        }
    }
}
