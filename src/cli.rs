use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "kira-proteoqc", version, about = "ProteoQC scaffolding CLI")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    Run(RunArgs),
    Geneset(GenesetArgs),
    Validate(ValidateArgs),
}

#[derive(Debug, Args)]
pub struct RunArgs {
    #[arg(long, num_args = 1.., help = "Input directory (10x MTX) or .h5ad file (repeatable)")]
    pub input: Vec<PathBuf>,

    #[arg(long)]
    pub out: PathBuf,

    #[arg(long, value_enum)]
    pub mode: ModeArg,

    #[arg(long, default_value_t = false)]
    pub timecourse: bool,

    #[arg(long)]
    pub geneset: Option<PathBuf>,

    #[arg(long, default_value_t = false)]
    pub no_log1p: bool,

    #[arg(long, default_value_t = false)]
    pub json: bool,

    #[arg(long, default_value_t = false)]
    pub tsv: bool,

    #[arg(long, default_value_t = 0, help = "Number of threads (0 = auto)")]
    pub threads: usize,

    #[arg(long, default_value_t = 4096, help = "Cache block size (cells)")]
    pub cache_block: usize,

    #[arg(
        long,
        default_value_t = false,
        help = "Enable prefetch (requires feature 'prefetch')"
    )]
    pub prefetch: bool,

    #[arg(
        long,
        default_value = "off",
        help = "Enable fusion mode: off|proteo|mito+proteo"
    )]
    pub fusion: String,

    #[arg(long, value_enum, default_value_t = RunModeArg::Standalone)]
    pub run_mode: RunModeArg,

    #[arg(long, help = "Path to shared cache file (kira-organelle.bin)")]
    pub cache: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct GenesetArgs {
    #[command(subcommand)]
    pub command: GenesetCommand,
}

#[derive(Debug, Subcommand)]
pub enum GenesetCommand {
    Show(GenesetShowArgs),
}

#[derive(Debug, Args)]
pub struct ValidateArgs {
    #[arg(long, help = "Input directory (10x MTX) or .h5ad file")]
    pub input: PathBuf,

    #[arg(long, value_enum, default_value_t = RunModeArg::Standalone)]
    pub run_mode: RunModeArg,
}

#[derive(Debug, Args)]
pub struct GenesetShowArgs {
    #[arg(long, help = "Optional input to resolve coverage")]
    pub input: Option<PathBuf>,

    #[arg(long, help = "Optional geneset TSV to overlay on built-in sets")]
    pub geneset: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ModeArg {
    Cell,
    Sample,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum RunModeArg {
    Standalone,
    Pipeline,
}
