//! Shared helpers for `pgr pl` pipeline subcommands.
//!
//! Pure pipeline orchestration logic (no clap dependency): workflow context,
//! path resolution, and external-tool driver functions (FastK / Profex / spanr).

mod ctx;
mod repeat;

pub use ctx::{abs_path_or_stdout, PipelineCtx};
pub use repeat::{run_profex_per_chr, run_repeat_pipeline, run_repeat_spanr_pipeline, RepeatOpts};
