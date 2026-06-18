use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
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
    min_control_non_g_frac: f64,
    /// Minimum total coverage required at this position in the RNA sample
    #[arg(long, default_value_t = 20)]
    min_rna_coverage: u32,
    /// Minimum edit (G) fraction required in the RNA sample
    #[arg(long, default_value_t = 0.05)]
    min_rna_edit_frac: f64,
}

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    let ctrlf = File::open(&cli.control_path)?;
    let rnaf = File::open(&cli.rna_path)?;
    let mut rnabuf = BufReader::new(&rnaf);
    let outfile = File::create(&cli.output_path)?;
    let mut writer = BufWriter::new(outfile);
    // Trust me bro
    //// SAFETY: ctrl file is written once by sam_to_map and not modified
    let ctrl_map = unsafe { Mmap::map(&ctrlf)? };
    let pos_records = std::iter::from_fn(move || {
        let mut buf = [0u8; 12];
        match rnabuf.read_exact(&mut buf) {
            Ok(()) => Some(buf),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(e) => panic!("read error {e}"),
        }
    });
    let pos_records = pos_records.map(|buf| {
        let pos = u32::from_le_bytes(buf[..4].try_into().unwrap());
        let g_count = u32::from_le_bytes(buf[4..8].try_into().unwrap());
        let other = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        (pos, g_count, other)
    });
    pos_records.for_each(|(pos, rna_g, rna_other)| {
        if let Some((ctrl_g, ctrl_other)) = lookup_control(&ctrl_map, pos) {
            let rna_tot = rna_g + rna_other;
            let ctrl_total = ctrl_g + ctrl_other;
            let rna_edit_frac = rna_g as f64 / rna_tot as f64;
            let ctrl_edit_frac = ctrl_g as f64 / ctrl_total as f64;
            let ctrl_nonedit_frac = ctrl_other as f64 / ctrl_total as f64;

            if ctrl_total > cli.min_control_coverage
                && ctrl_edit_frac < cli.max_control_edit_frac
                && ctrl_nonedit_frac >= cli.min_control_non_g_frac
                && rna_tot >= cli.min_rna_coverage
                && rna_edit_frac >= cli.min_rna_edit_frac
            {
                writer.write_all(&pos.to_le_bytes()).unwrap();
                writer.write_all(&rna_g.to_le_bytes()).unwrap();
                writer.write_all(&ctrl_g.to_le_bytes()).unwrap();
            }
        }
    });
    writer.flush()?;
    Ok(())
}

fn lookup_control(ctrl_mmap: &[u8], target: u32) -> Option<(u32, u32)> {
    let n_records = ctrl_mmap.len() / 12;
    if n_records == 0 {
        return None;
    }

    let mut low = 0;
    let mut high = n_records - 1;

    while low <= high {
        let mid = (low + high) / 2;
        let offset = mid * 12;

        let pos = u32::from_le_bytes(ctrl_mmap[offset..offset + 4].try_into().unwrap());

        if pos == target {
            let g = u32::from_le_bytes(ctrl_mmap[offset + 4..offset + 8].try_into().unwrap());
            let other = u32::from_le_bytes(ctrl_mmap[offset + 8..offset + 12].try_into().unwrap());
            return Some((g, other));
        } else if target < pos {
            if mid == 0 {
                return None;
            }
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }
    None
}
