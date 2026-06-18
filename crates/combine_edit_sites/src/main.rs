use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::path::PathBuf;

#[derive(Parser)]
struct Cli {
    /// Input binary files (one per replicate)
    inputs: Vec<PathBuf>,
    /// Output binary file
    #[arg(long)]
    output: PathBuf,
    /// Combination mode: and/or
    #[arg(long)]
    mode: Mode,
}

#[derive(clap::ValueEnum, Clone)]
enum Mode {
    And,
    Or,
}

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    let files: Vec<File> = cli
        .inputs
        .iter()
        .map(|path| File::open(&path).expect("Cannot open input file"))
        .collect();
    // SAFETY: input files is written once by sam_to_map and not modified
    let mmaps: Vec<Mmap> = files
        .iter()
        .map(|f| unsafe { Mmap::map(f).expect("couldn't mmap file") })
        .collect();
    Ok(())
}

fn current_pos(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * 12;
    if offset + 4 > mmap.len() {
        None // exhausted
    } else {
        Some(u32::from_le_bytes(
            mmap[offset..offset + 4].try_into().unwrap(),
        ))
    }
}

