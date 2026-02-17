use clap::Parser;
use kira_proteoqc::cli::{Cli, Commands, RunModeArg};

#[test]
fn run_mode_defaults_to_standalone() {
    let cli = Cli::parse_from([
        "kira-proteoqc",
        "run",
        "--input",
        "data",
        "--out",
        "out",
        "--mode",
        "cell",
    ]);
    match cli.command {
        Commands::Run(args) => assert_eq!(args.run_mode, RunModeArg::Standalone),
        _ => panic!("expected run command"),
    }
}

#[test]
fn run_mode_pipeline_is_accepted() {
    let cli = Cli::parse_from([
        "kira-proteoqc",
        "run",
        "--input",
        "data",
        "--out",
        "out",
        "--mode",
        "cell",
        "--run-mode",
        "pipeline",
    ]);
    match cli.command {
        Commands::Run(args) => assert_eq!(args.run_mode, RunModeArg::Pipeline),
        _ => panic!("expected run command"),
    }
}
