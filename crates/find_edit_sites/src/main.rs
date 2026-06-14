use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::Read;
use std::thread::panicking;
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
    let ctrlf = File::open(&cli.control_path)?;
    let rnaf = File::open(&cli.rna_path)?;
    let mut rnabuf = BufReader::new(&rnaf);
    // Trust me bro
    let ctrl_map = unsafe { Mmap::map(&ctrlf)? };
    let pos_records = std::iter::from_fn(move || {
        let mut buf = [0u8; 16];
        match rnabuf.read_exact(&mut buf) {
            Ok(()) => Some(buf),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(e) => panic!("read error {e}"),
        }
    });
    let pos_records = pos_records.map(|buf| {
        let pos = u32::from_le_bytes(buf[..4].try_into().unwrap());
        let a_count = u32::from_le_bytes(buf[4..8].try_into().unwrap());
        let g_count = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        let other = u32::from_le_bytes(buf[12..16].try_into().unwrap());
        (pos, a_count, g_count, other)
    });

    Ok(())
}

fn get_position(slice: &[u8]) -> u32 {
    let bytes: [u8; 4] = slice.try_into().unwrap();
    u32::from_le_bytes(bytes)
}

fn lookup_control(ctrl_mmap: &[u8], target: u32) -> Option<(u32, u32, u32)> {
    let n_records = ctrl_mmap.len() / 16;
    if n_records == 0 {
        return None;
    }

    let mut low = 0;
    let mut high = n_records - 1;

    while low <= high {
        let mid = (low + high) / 2;
        let offset = mid * 16;

        let pos = u32::from_le_bytes(ctrl_mmap[offset..offset + 4].try_into().unwrap());

        if pos == target {
            let a = u32::from_le_bytes(ctrl_mmap[offset + 4..offset + 8].try_into().unwrap());
            let g = u32::from_le_bytes(ctrl_mmap[offset + 8..offset + 12].try_into().unwrap());
            let other = u32::from_le_bytes(ctrl_mmap[offset + 12..offset + 16].try_into().unwrap());
            return Some((a, g, other));
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
