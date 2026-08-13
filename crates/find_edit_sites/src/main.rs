use clap::Parser;
use memmap2::Mmap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::{io::BufReader, path::PathBuf};

#[derive(Parser)]
struct Cli {
    rna_path: PathBuf,
    control_path: PathBuf,
    output_path: PathBuf,
    /// Sorted synthetic gene-span GTF (PREPARE_GENE_SPANS output) — used to
    /// resolve which strand/channel (A->G vs T->C) is the real edit signal
    /// at each position, and to restrict candidates to positions inside a gene
    gene_spans_path: PathBuf,
    chr_list_path: PathBuf,
    #[arg(long, default_value_t = 10)]
    min_control_coverage: u32,
    /// Maximum allowed edit fraction in the control sample (filters out SNPs)
    #[arg(long, default_value_t = 0.005)]
    max_control_edit_frac: f64,
    #[arg(long, default_value_t = 0.8)]
    min_control_non_g_frac: f64,
    /// Minimum edit-signal read count required in the RNA sample (legacy: > 0)
    #[arg(long, default_value_t = 1)]
    min_rna_edit_count: u32,
}

#[derive(Clone, Copy)]
struct GeneSpan {
    start: u32,
    end: u32,
    strand_plus: bool,
}

// 21-byte sam_to_map record: chr_id(1) pos(4) a_to_g(4) a_match(4) t_to_c(4) t_match(4)
const REC_SIZE: usize = 21;

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();

    let chr_names: Vec<String> = std::fs::read_to_string(&cli.chr_list_path)
        .unwrap()
        .lines()
        .map(|s| s.to_string())
        .collect();
    let chr_lookup: HashMap<&str, u8> = chr_names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i as u8))
        .collect();

    // group gene spans per chromosome, sorted by start, for lookup
    let mut gene_spans: HashMap<u8, Vec<GeneSpan>> = HashMap::new();
    for line in std::fs::read_to_string(&cli.gene_spans_path)
        .unwrap()
        .lines()
    {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            continue;
        }
        let Some(&chr_id) = chr_lookup.get(f[0]) else {
            continue;
        };
        let start: u32 = f[3].parse().unwrap();
        let end: u32 = f[4].parse().unwrap();
        let strand_plus = f[6] == "+";
        gene_spans.entry(chr_id).or_default().push(GeneSpan {
            start,
            end,
            strand_plus,
        });
    }
    for spans in gene_spans.values_mut() {
        spans.sort_by_key(|s| s.start);
    }

    let ctrlf = File::open(&cli.control_path)?;
    let rnaf = File::open(&cli.rna_path)?;
    let mut rnabuf = BufReader::new(&rnaf);
    let outfile = File::create(&cli.output_path)?;
    let mut writer = BufWriter::new(outfile);

    // SAFETY: ctrl file is written once by sam_to_map and not modified
    let ctrl_map = unsafe { Mmap::map(&ctrlf)? };

    let pos_records = std::iter::from_fn(move || {
        let mut buf = [0u8; REC_SIZE];
        match rnabuf.read_exact(&mut buf) {
            Ok(()) => Some(buf),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(e) => panic!("read error {e}"),
        }
    });
    let pos_records = pos_records.map(|buf| {
        let chr_id = buf[0];
        let pos = u32::from_le_bytes(buf[1..5].try_into().unwrap());
        let a_to_g = u32::from_le_bytes(buf[5..9].try_into().unwrap());
        let a_match = u32::from_le_bytes(buf[9..13].try_into().unwrap());
        let t_to_c = u32::from_le_bytes(buf[13..17].try_into().unwrap());
        let t_match = u32::from_le_bytes(buf[17..21].try_into().unwrap());
        (chr_id, pos, a_to_g, a_match, t_to_c, t_match)
    });

    pos_records.for_each(
        |(chr_id, pos, rna_a_to_g, rna_a_match, rna_t_to_c, rna_t_match)| {
            let Some(spans) = gene_spans.get(&chr_id) else {
                return;
            };

            // interval-stabbing: find every gene span overlapping this position.
            // Linear scan restricted to this chromosome's genes — fine at
            // test-genome scale, but revisit with an interval tree if this
            // becomes a bottleneck on a full genome.
            for gene in spans.iter().filter(|g| g.start <= pos && pos <= g.end) {
                let (rna_g, rna_other) = if gene.strand_plus {
                    (rna_a_to_g, rna_a_match)
                } else {
                    (rna_t_to_c, rna_t_match)
                };
                if rna_g < cli.min_rna_edit_count {
                    continue;
                }

                let Some((ctrl_g, ctrl_other)) =
                    lookup_control(&ctrl_map, chr_id, pos, gene.strand_plus)
                else {
                    continue;
                };
                let rna_tot = rna_g + rna_other;
                let ctrl_total = ctrl_g + ctrl_other;
                if rna_tot == 0 || ctrl_total == 0 {
                    continue;
                }
                let ctrl_edit_frac = ctrl_g as f64 / ctrl_total as f64;
                let ctrl_nonedit_frac = ctrl_other as f64 / ctrl_total as f64;

                let passes = ctrl_total >= cli.min_control_coverage
                    && ctrl_edit_frac < cli.max_control_edit_frac
                    && ctrl_nonedit_frac >= cli.min_control_non_g_frac;

                if passes {
                    writer.write_all(&[chr_id]).unwrap();
                    writer.write_all(&pos.to_le_bytes()).unwrap();
                    writer.write_all(&rna_g.to_le_bytes()).unwrap();
                    writer.write_all(&rna_other.to_le_bytes()).unwrap();
                    writer.write_all(&ctrl_g.to_le_bytes()).unwrap();
                    writer.write_all(&ctrl_other.to_le_bytes()).unwrap();
                }
            }
        },
    );

    writer.flush()?;
    Ok(())
}

fn lookup_control(
    ctrl_mmap: &[u8],
    target_chr: u8,
    target_pos: u32,
    strand_plus: bool,
) -> Option<(u32, u32)> {
    let n_records = ctrl_mmap.len() / REC_SIZE;
    if n_records == 0 {
        return None;
    }
    let mut low = 0;
    let mut high = n_records - 1;
    while low <= high {
        let mid = (low + high) / 2;
        let offset = mid * REC_SIZE;
        let chr_id = ctrl_mmap[offset];
        let pos = u32::from_le_bytes(ctrl_mmap[offset + 1..offset + 5].try_into().unwrap());
        match (chr_id, pos).cmp(&(target_chr, target_pos)) {
            std::cmp::Ordering::Equal => {
                let (g, other) = if strand_plus {
                    let g =
                        u32::from_le_bytes(ctrl_mmap[offset + 5..offset + 9].try_into().unwrap());
                    let other =
                        u32::from_le_bytes(ctrl_mmap[offset + 9..offset + 13].try_into().unwrap());
                    (g, other)
                } else {
                    let g =
                        u32::from_le_bytes(ctrl_mmap[offset + 13..offset + 17].try_into().unwrap());
                    let other =
                        u32::from_le_bytes(ctrl_mmap[offset + 17..offset + 21].try_into().unwrap());
                    (g, other)
                };
                return Some((g, other));
            }
            std::cmp::Ordering::Greater => {
                if mid == 0 {
                    return None;
                }
                high = mid - 1;
            }
            std::cmp::Ordering::Less => low = mid + 1,
        }
    }
    None
}
