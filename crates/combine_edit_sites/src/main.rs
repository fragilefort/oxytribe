use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
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
    let outfile = File::create(&cli.output)?;
    let mut writer = BufWriter::new(outfile);
    let mmaps: Vec<Mmap> = files
        .iter()
        .map(|f| unsafe { Mmap::map(f).expect("couldn't mmap file") })
        .collect();
    let pointers = vec![0; mmaps.len()];
    let n_files = mmaps.len();
    while pointers
        .iter()
        .enumerate()
        .any(|(i, &p)| current_pos(&mmaps[i], p).is_some())
    {
        let min_pos = pointers
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| current_pos(&mmaps[i], p))
            .min()
            .unwrap();

        let contributors: Vec<usize> = pointers
            .iter()
            .enumerate()
            .filter(|&(i, &p)| current_pos(&mmaps[i], p) == Some(min_pos))
            .map(|(i, _)| i)
            .collect();

        let emit = match cli.mode {
            Mode::And => contributors.len() == n_files,
            Mode::Or => true,
        };
        if emit {
            // compute mean g and other across contributors
            let mean_g = contributors
                .iter()
                .map(|&i| get_g(&mmaps[i], pointers[i]))
                .sum::<u32>()
                / contributors.len() as u32;
            let mean_other = contributors
                .iter()
                .map(|&i| get_other(&mmaps[i], pointers[i]))
                .sum::<u32>()
                / contributors.len() as u32;

            writer.write_all(&min_pos.to_le_bytes()).unwrap();
            writer.write_all(&mean_g.to_le_bytes()).unwrap();
            writer.write_all(&mean_other.to_le_bytes()).unwrap();
        }
        for i in &contributors {
            pointers[*i] += 1;
        }
    }
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

fn get_g(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * 12;
    if offset + 12 > mmap.len() {
        None
    } else {
        Some(u32::from_le_bytes(
            mmap[offset + 4..offset + 8].try_into().unwrap(),
        ))
    }
}

fn get_other(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * 12;
    if offset + 12 > mmap.len() {
        None
    } else {
        Some(u32::from_le_bytes(
            mmap[offset + 8..offset + 12].try_into().unwrap(),
        ))
    }
}
