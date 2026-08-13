use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser)]
struct Cli {
    inputs: Vec<PathBuf>,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value = "or")]
    mode: Mode,
    /// Chromosome list (required when outputting TSV)
    #[arg(long)]
    chr_list: Option<PathBuf>,
    /// Pass if pooling raw 13-byte control maps (outputs a 13-byte .bin)
    #[arg(long)]
    raw: bool,
}

#[derive(clap::ValueEnum, Clone)]
enum Mode {
    And,
    Or,
}

// Controls: chr_id(1) pos(4) g(4) other(4) = 13 bytes
// Filtered treatments: chr_id(1) pos(4) strand(1) rna_g(4) rna_other(4) ctrl_g(4) ctrl_other(4) = 22 bytes
const RAW_REC_SIZE: usize = 13;
const FILTERED_REC_SIZE: usize = 22;

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();

    let record_size: usize = if cli.raw {
        RAW_REC_SIZE
    } else {
        FILTERED_REC_SIZE
    };

    let chr_names: Vec<String> = cli
        .chr_list
        .map(|p| {
            std::fs::read_to_string(p)
                .unwrap()
                .lines()
                .map(String::from)
                .collect()
        })
        .unwrap_or_default();

    let files: Vec<File> = cli
        .inputs
        .iter()
        .map(|path| File::open(path).expect("Cannot open input file"))
        .collect();

    let outfile = File::create(&cli.output)?;
    let mut writer = BufWriter::new(outfile);

    if !cli.raw {
        writeln!(
            writer,
            "#chr\tstart\tend\trna_G\trna_total\tctrl_G\tctrl_total\trna_edit_frac"
        )?;
    }

    let mmaps: Vec<Mmap> = files
        .iter()
        .map(|f| unsafe { Mmap::map(f).expect("couldn't mmap file") })
        .collect();

    let mut pointers: Vec<usize> = vec![0; mmaps.len()];
    let n_files = mmaps.len();

    // key = (chr_id, pos, strand). strand is always 0 in raw mode — controls
    // have no strand field, and pooling them isn't strand-aware, so treating
    // strand as a constant there is a correct no-op, not a workaround.
    let current_key = |mmap: &[u8], ptr: usize| -> Option<(u8, u32, u8)> {
        let offset = ptr * record_size;
        if offset + record_size > mmap.len() {
            None
        } else {
            let chr_id = mmap[offset];
            let pos = u32::from_le_bytes(mmap[offset + 1..offset + 5].try_into().unwrap());
            let strand = if cli.raw { 0 } else { mmap[offset + 5] };
            Some((chr_id, pos, strand))
        }
    };

    while pointers
        .iter()
        .enumerate()
        .any(|(i, &p)| current_key(&mmaps[i], p).is_some())
    {
        let min_key = pointers
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| current_key(&mmaps[i], p))
            .min()
            .unwrap();

        let contributors: Vec<usize> = pointers
            .iter()
            .enumerate()
            .filter(|&(i, &p)| current_key(&mmaps[i], p) == Some(min_key))
            .map(|(i, _)| i)
            .collect();

        let emit = match cli.mode {
            Mode::And => contributors.len() == n_files,
            Mode::Or => true,
        };

        if emit {
            let (g_off, other_off) = if cli.raw { (5, 9) } else { (6, 10) };

            let sum_g = contributors
                .iter()
                .map(|&i| {
                    let off = pointers[i] * record_size;
                    u32::from_le_bytes(mmaps[i][off + g_off..off + g_off + 4].try_into().unwrap())
                })
                .sum::<u32>();

            let sum_other = contributors
                .iter()
                .map(|&i| {
                    let off = pointers[i] * record_size;
                    u32::from_le_bytes(
                        mmaps[i][off + other_off..off + other_off + 4]
                            .try_into()
                            .unwrap(),
                    )
                })
                .sum::<u32>();

            let sum_total = sum_g + sum_other;

            if cli.raw {
                if sum_total > 0 {
                    let (chr_id, pos, _strand) = min_key;
                    writer.write_all(&[chr_id])?;
                    writer.write_all(&pos.to_le_bytes())?;
                    writer.write_all(&sum_g.to_le_bytes())?;
                    writer.write_all(&sum_other.to_le_bytes())?;
                }
            } else if sum_total > 0 {
                let first = contributors[0];
                let off = pointers[first] * record_size;
                let ctrl_g =
                    u32::from_le_bytes(mmaps[first][off + 14..off + 18].try_into().unwrap());
                let ctrl_other =
                    u32::from_le_bytes(mmaps[first][off + 18..off + 22].try_into().unwrap());
                let ctrl_total = ctrl_g + ctrl_other;

                let (chr_id, pos, _strand) = min_key;
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
                )?;
            }
        }

        for i in &contributors {
            pointers[*i] += 1;
        }
    }

    writer.flush()?;
    Ok(())
}
