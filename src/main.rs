use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;
use tracing_subscriber::EnvFilter;

use kira_proteoqc::cli::{Cli, Commands, ModeArg, RunModeArg};
use kira_proteoqc::ctx::{Ctx, RunMode};
use kira_proteoqc::geneset;
use kira_proteoqc::io;
use kira_proteoqc::pipeline::Pipeline;
use kira_proteoqc::pipeline::stage0_scaffold::Stage0Scaffold;
use kira_proteoqc::pipeline::stage1_input::Stage1Input;
use kira_proteoqc::pipeline::stage2_h5ad::Stage2H5ad;
use kira_proteoqc::pipeline::stage3_expr_cache::Stage3ExprCache;
use kira_proteoqc::pipeline::stage4_geneset::Stage4Geneset;
use kira_proteoqc::pipeline::stage5_math::Stage5Math;
use kira_proteoqc::pipeline::stage6_axes::Stage6Axes;
use kira_proteoqc::pipeline::stage7_integrate::Stage7Integrate;
use kira_proteoqc::pipeline::stage8_risk::Stage8Risk;
use kira_proteoqc::pipeline::stage9_timecourse::Stage9Timecourse;
use kira_proteoqc::pipeline::stage10_output::Stage10Output;
use kira_proteoqc::schema::v1::Mode;
use kira_proteoqc::scores::TimepointSummary;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();
    tracing::info!(
        simd_backend = %kira_proteoqc::simd::backend_name(),
        "simd backend selected"
    );

    match cli.command {
        Commands::Run(args) => {
            let mode = match args.mode {
                ModeArg::Cell => Mode::Cell,
                ModeArg::Sample => Mode::Sample,
            };
            let log1p = !args.no_log1p;

            if args.timecourse && args.input.len() < 2 {
                anyhow::bail!("--timecourse requires at least 2 --input values");
            }
            if !args.timecourse && args.input.len() != 1 {
                anyhow::bail!("multiple --input requires --timecourse");
            }

            if args.timecourse {
                let ordered_inputs = order_timecourse_inputs(&args.input);
                let mut master_ctx = Ctx::new(
                    ordered_inputs[0].clone(),
                    args.out.clone(),
                    mode.clone(),
                    args.timecourse,
                    args.geneset.clone(),
                    log1p,
                    args.json,
                    args.tsv,
                    env!("CARGO_PKG_VERSION"),
                );
                master_ctx.threads = args.threads;
                master_ctx.cache_block = args.cache_block;
                master_ctx.prefetch = args.prefetch;
                master_ctx.fusion = args.fusion.clone();
                master_ctx.run_mode = match args.run_mode {
                    RunModeArg::Standalone => RunMode::Standalone,
                    RunModeArg::Pipeline => RunMode::Pipeline,
                };
                master_ctx.cache_override = args.cache.clone();

                for input in ordered_inputs {
                    let label = label_from_path(&input);
                    let out_dir = master_ctx.output.out_dir.join(&label);
                    let mut ctx = Ctx::new(
                        input,
                        out_dir,
                        mode.clone(),
                        args.timecourse,
                        args.geneset.clone(),
                        log1p,
                        args.json,
                        args.tsv,
                        env!("CARGO_PKG_VERSION"),
                    );
                    ctx.threads = args.threads;
                    ctx.cache_block = args.cache_block;
                    ctx.prefetch = args.prefetch;
                    ctx.fusion = args.fusion.clone();
                    ctx.run_mode = match args.run_mode {
                        RunModeArg::Standalone => RunMode::Standalone,
                        RunModeArg::Pipeline => RunMode::Pipeline,
                    };
                    ctx.cache_override = args.cache.clone();
                    let pipeline = Pipeline::new(vec![
                        Box::new(Stage0Scaffold::new()),
                        Box::new(Stage1Input::new()),
                        Box::new(Stage2H5ad::new()),
                        Box::new(Stage3ExprCache::new()),
                        Box::new(Stage4Geneset::new()),
                        Box::new(Stage5Math::new()),
                        Box::new(Stage6Axes::new()),
                        Box::new(Stage7Integrate::new()),
                        Box::new(Stage8Risk::new()),
                        Box::new(Stage10Output::new()),
                    ]);
                    pipeline.run(&mut ctx)?;
                    let tp = build_timepoint_summary(&ctx, label)?;
                    master_ctx.timecourse_points.push(tp);
                }

                let pipeline = Pipeline::new(vec![Box::new(Stage9Timecourse::new())]);
                pipeline.run(&mut master_ctx)?;
                print_timecourse_summary(&master_ctx);
            } else {
                let input = args.input.into_iter().next().unwrap();
                let mut ctx = Ctx::new(
                    input,
                    args.out,
                    mode,
                    args.timecourse,
                    args.geneset,
                    log1p,
                    args.json,
                    args.tsv,
                    env!("CARGO_PKG_VERSION"),
                );
                ctx.threads = args.threads;
                ctx.cache_block = args.cache_block;
                ctx.prefetch = args.prefetch;
                ctx.fusion = args.fusion;
                ctx.run_mode = match args.run_mode {
                    RunModeArg::Standalone => RunMode::Standalone,
                    RunModeArg::Pipeline => RunMode::Pipeline,
                };
                ctx.cache_override = args.cache.clone();

                let pipeline = Pipeline::new(vec![
                    Box::new(Stage0Scaffold::new()),
                    Box::new(Stage1Input::new()),
                    Box::new(Stage2H5ad::new()),
                    Box::new(Stage3ExprCache::new()),
                    Box::new(Stage4Geneset::new()),
                    Box::new(Stage5Math::new()),
                    Box::new(Stage6Axes::new()),
                    Box::new(Stage7Integrate::new()),
                    Box::new(Stage8Risk::new()),
                    Box::new(Stage9Timecourse::new()),
                    Box::new(Stage10Output::new()),
                ]);
                pipeline.run(&mut ctx)?;

                print_summary(&ctx)?;
            }
        }
        Commands::Geneset(args) => match args.command {
            kira_proteoqc::cli::GenesetCommand::Show(show) => {
                handle_geneset_show(show)?;
            }
        },
        Commands::Validate(args) => {
            let mut ctx = Ctx::new(
                args.input,
                PathBuf::from("."),
                Mode::Cell,
                false,
                None,
                true,
                false,
                false,
                env!("CARGO_PKG_VERSION"),
            );
            ctx.run_mode = match args.run_mode {
                RunModeArg::Standalone => RunMode::Standalone,
                RunModeArg::Pipeline => RunMode::Pipeline,
            };

            let pipeline = Pipeline::new(vec![
                Box::new(Stage1Input::new()),
                Box::new(Stage2H5ad::new()),
                Box::new(Stage3ExprCache::new()),
                Box::new(Stage4Geneset::new()),
                Box::new(Stage5Math::new()),
                Box::new(Stage6Axes::new()),
                Box::new(Stage7Integrate::new()),
                Box::new(Stage8Risk::new()),
                Box::new(Stage9Timecourse::new()),
            ]);
            pipeline.run(&mut ctx)?;

            print_validate_summary(&ctx);
        }
    }

    Ok(())
}

fn print_summary(ctx: &Ctx) -> Result<()> {
    let summary = io::summary::format_summary(ctx)?;
    print!("{}", summary);
    if !ctx.warnings.is_empty() {
        println!("warnings:");
        for warning in &ctx.warnings {
            println!("- {}", warning);
        }
    }
    Ok(())
}

fn print_validate_summary(ctx: &Ctx) {
    println!("kira-proteoqc validate ok");
    println!("genes: {}", ctx.genes.len());
    println!("cells: {}", ctx.cells.len());
    println!("nnz: {}", ctx.nnz);
    if !ctx.warnings.is_empty() {
        println!("warnings:");
        for warning in &ctx.warnings {
            println!("- {}", warning);
        }
    }
}

fn handle_geneset_show(args: kira_proteoqc::cli::GenesetShowArgs) -> Result<()> {
    if let Some(input) = args.input {
        let mut ctx = Ctx::new(
            input,
            PathBuf::from("."),
            Mode::Cell,
            false,
            args.geneset,
            true,
            false,
            false,
            env!("CARGO_PKG_VERSION"),
        );
        let pipeline = Pipeline::new(vec![
            Box::new(Stage1Input::new()),
            Box::new(Stage2H5ad::new()),
            Box::new(Stage4Geneset::new()),
        ]);
        pipeline.run(&mut ctx)?;
        print_geneset_with_coverage(&ctx)?;
        return Ok(());
    }

    let mut collection = geneset::load_builtin()?;
    if let Some(path) = args.geneset {
        let user_defs = geneset::load_user(&path)?;
        collection.defs = geneset::merge_defs(collection.defs, user_defs);
    }
    print_geneset_list(&collection);
    Ok(())
}

fn print_geneset_list(collection: &geneset::GenesetCollection) {
    println!("genesets (version {}):", collection.version);
    for def in &collection.defs {
        println!("{}\t{}\t{}", def.id, def.axis, def.genes.len());
    }
}

fn print_geneset_with_coverage(ctx: &Ctx) -> Result<()> {
    let collection = ctx.genesets.as_ref().context("genesets not resolved")?;
    println!("genesets (version {}):", collection.version);
    for gs in &collection.resolved {
        let found = gs.gene_ids.len();
        let total = gs.total;
        let fraction = if total == 0 {
            0.0
        } else {
            found as f64 / total as f64
        };
        println!(
            "{}\t{}\t{}\t{}/{}\t{:.4}",
            gs.id, gs.axis, total, found, total, fraction
        );
    }
    Ok(())
}

fn build_timepoint_summary(ctx: &Ctx, label: String) -> Result<TimepointSummary> {
    let axis = ctx
        .axis_raw
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("axis raw missing"))?;
    let integrated = ctx
        .integrated_scores
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("integrated scores missing"))?;

    let pcs = mean_vec(&axis.pcs)?;
    let cls = mean_vec(&axis.cls)?;
    let utp = mean_vec(&axis.utp)?;
    let pfs = mean_vec(&integrated.pfs_raw)?;
    let pii = mean_vec(&integrated.pii_raw)?;

    Ok(TimepointSummary {
        label,
        pfs,
        pii,
        pcs,
        cls,
        utp,
    })
}

fn mean_vec(values: &[f32]) -> Result<f32> {
    if values.is_empty() {
        return Ok(0.0);
    }
    let mut sum = 0.0;
    for v in values {
        if v.is_nan() {
            anyhow::bail!("NaN encountered in timecourse summary");
        }
        sum += *v;
    }
    Ok(sum / values.len() as f32)
}

fn order_timecourse_inputs(inputs: &[PathBuf]) -> Vec<PathBuf> {
    let mut labels: Vec<(String, PathBuf)> = inputs
        .iter()
        .map(|p| (label_from_path(p), p.clone()))
        .collect();
    let all_have_t = labels.iter().all(|(l, _)| has_time_token(l));
    if all_have_t {
        labels.sort_by(|a, b| a.0.cmp(&b.0));
    }
    labels.into_iter().map(|(_, p)| p).collect()
}

fn has_time_token(label: &str) -> bool {
    if let Some(pos) = label.find("_T") {
        let rest = &label[pos + 2..];
        return !rest.is_empty() && rest.chars().take_while(|c| c.is_ascii_digit()).count() > 0;
    }
    false
}

fn label_from_path(path: &PathBuf) -> String {
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("input")
        .to_string()
}

fn print_timecourse_summary(ctx: &Ctx) {
    if let Some(tc) = &ctx.timecourse_result {
        println!("timecourse trajectory: {}", tc.trajectory);
        for d in &tc.deltas {
            println!(
                "{} -> {}: dPFS={:.4} dPII={:.4} dPCS={:.4} dCLS={:.4}",
                d.from, d.to, d.delta_pfs, d.delta_pii, d.delta_pcs, d.delta_cls
            );
        }
    }
}
