use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
struct Cli {
    /// Input binary files (one per replicate), each already filtered by find_edit_sites
    inputs: Vec<PathBuf>,
    /// Output TSV file
    #[arg(long)]
    output: PathBuf,
    /// Combination mode: and/or
    #[arg(long)]
    mode: Mode,
    /// Chromosome name list, same file used by find_edit_sites, to resolve chr_id -> chr name
    #[arg(long)]
    chr_list: PathBuf,
}

#[derive(clap::ValueEnum, Clone)]
enum Mode {
    And,
    Or,
}

const RECORD_SIZE: usize = 21; // chr_id(1) + pos(4) + rna_g(4) + rna_other(4) + ctrl_g(4) + ctrl_other(4)

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    let chr_names: Vec<String> = std::fs::read_to_string(&cli.chr_list)
        .unwrap()
        .lines()
        .map(|s| s.to_string())
        .collect();

    let files: Vec<File> = cli
        .inputs
        .iter()
        .map(|path| File::open(path).expect("Cannot open input file"))
        .collect();

    let outfile = File::create(&cli.output)?;
    let mut writer = BufWriter::new(outfile);
    writeln!(
        writer,
        "#chr\tstart\tend\trna_G\trna_total\tctrl_G\tctrl_total\trna_edit_frac"
    )?;

    // SAFETY: input files are written once by find_edit_sites and not modified
    let mmaps: Vec<Mmap> = files
        .iter()
        .map(|f| unsafe { Mmap::map(f).expect("couldn't mmap file") })
        .collect();

    let mut pointers: Vec<usize> = vec![0; mmaps.len()];
    let n_files = mmaps.len();

    while pointers
        .iter()
        .enumerate()
        .any(|(i, &p)| current_pos(&mmaps[i], p).is_some())
    {
        let min_key = pointers
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| current_pos(&mmaps[i], p))
            .min()
            .unwrap();

        let contributors: Vec<usize> = pointers
            .iter()
            .enumerate()
            .filter(|&(i, &p)| current_pos(&mmaps[i], p) == Some(min_key))
            .map(|(i, _)| i)
            .collect();

        let emit = match cli.mode {
            Mode::And => contributors.len() == n_files,
            Mode::Or => true,
        };

        if emit {
            let sum_g = contributors
                .iter()
                .map(|&i| get_rna_g(&mmaps[i], pointers[i]).unwrap())
                .sum::<u32>();
            let sum_other = contributors
                .iter()
                .map(|&i| get_rna_other(&mmaps[i], pointers[i]).unwrap())
                .sum::<u32>();
            let sum_total = sum_g + sum_other;

            // control counts are the same underlying control sample regardless
            // of which replicate file we read them from — take the first
            // contributor's values rather than summing.
            let first = contributors[0];
            let ctrl_g = get_ctrl_g(&mmaps[first], pointers[first]).unwrap();
            let ctrl_other = get_ctrl_other(&mmaps[first], pointers[first]).unwrap();
            let ctrl_total = ctrl_g + ctrl_other;

            if sum_total > 0 {
                let (chr_id, pos) = min_key;
                let chr_name = &chr_names[chr_id as usize];
                let rna_edit_frac = sum_g as f64 / sum_total as f64;

                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                    chr_name,
                    pos - 1,
                    pos,
                    sum_g,
                    sum_total,
                    ctrl_g,
                    ctrl_total,
                    rna_edit_frac
                )
                .unwrap();
            }
        }

        for i in &contributors {
            pointers[*i] += 1;
        }
    }

    writer.flush()?;
    Ok(())
}

fn current_pos(mmap: &[u8], pointer: usize) -> Option<(u8, u32)> {
    let offset = pointer * RECORD_SIZE;
    if offset + RECORD_SIZE > mmap.len() {
        None
    } else {
        let chr_id = mmap[offset];
        let pos = u32::from_le_bytes(mmap[offset + 1..offset + 5].try_into().unwrap());
        Some((chr_id, pos))
    }
}

fn get_rna_g(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * RECORD_SIZE;
    if offset + RECORD_SIZE > mmap.len() {
        None
    } else {
        Some(u32::from_le_bytes(
            mmap[offset + 5..offset + 9].try_into().unwrap(),
        ))
    }
}

fn get_rna_other(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * RECORD_SIZE;
    if offset + RECORD_SIZE > mmap.len() {
        None
    } else {
        Some(u32::from_le_bytes(
            mmap[offset + 9..offset + 13].try_into().unwrap(),
        ))
    }
}

fn get_ctrl_g(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * RECORD_SIZE;
    if offset + RECORD_SIZE > mmap.len() {
        None
    } else {
        Some(u32::from_le_bytes(
            mmap[offset + 13..offset + 17].try_into().unwrap(),
        ))
    }
}

fn get_ctrl_other(mmap: &[u8], pointer: usize) -> Option<u32> {
    let offset = pointer * RECORD_SIZE;
    if offset + RECORD_SIZE > mmap.len() {
        None
    } else {
        Some(u32::from_le_bytes(
            mmap[offset + 17..offset + 21].try_into().unwrap(),
        ))
    }
}
