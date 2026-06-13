use clap::Parser;
use std::{io::BufReader, path::PathBuf};

#[derive(Parser)]
struct Cli {
    rna_path: PathBuf,
    control_path: PathBuf,
    output_path: PathBuf,
    #[arg(long, default_value_t = 9)]
    min_control_coverage: u32,
    /// Maximum allowed edit (G) fraction in the control sample (filters out SNPs)
    #[arg(long, default_value_t = 0.005)]
    max_control_edit_frac: f64,
    #[arg(long, default_value_t = 0.8)]
    min_control_nonedit_frac: f64,
    /// Minimum total coverage required at this position in the RNA sample
    #[arg(long, default_value_t = 20)]
    min_rna_coverage: u32,
    /// Minimum edit (G) fraction required in the RNA sample
    #[arg(long, default_value_t = 0.05)]
    min_rna_edit_frac: f64,
}
fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    let ctrlf = open(&cli.control_path)?;
    let rnaf = open(&cli.rna_path)?;

    let rnabuf = BufReader::new(rnaf);
    Ok(())
}
