use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::PathBuf;

static CIGAR_TABLE: [(u8, u8); 128] = {
    let mut table = [(0u8, 0u8); 128];
    table[b'M' as usize] = (1, 1);
    table[b'=' as usize] = (1, 1);
    table[b'X' as usize] = (1, 1);
    table[b'I' as usize] = (1, 0);
    table[b'S' as usize] = (1, 0);
    table[b'D' as usize] = (0, 1);
    table[b'N' as usize] = (0, 1);
    table[b'H' as usize] = (0, 0);
    table[b'P' as usize] = (0, 0);
    table
};

#[derive(Parser)]
struct Cli {
    input_sam: PathBuf,
    output_path: PathBuf,
    chr_list_path: PathBuf,
}

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

    let mut map: HashMap<(u8, u32), (u32, u32)> = HashMap::with_capacity(1_000_000);
    let f = File::open(&cli.input_sam).expect("Couldn't open sam file");
    let reader = BufReader::new(f);
    reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|line| !line.starts_with("@"))
        .for_each(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let chr = fields[2];
            let Some(&chr_id) = chr_lookup.get(&chr) else {
                return;
            };
            let start = fields[3].parse::<u32>().unwrap();
            let cigar = fields[5];
            let sequence = fields[9].as_bytes();
            let md_tag = extract_md_tag(&fields).unwrap();
            let cigar_per_base = parse_cigar(cigar)
                .flat_map(|(len, r, q)| std::iter::repeat((r, q)).take(len as usize));
            let mut ref_pos = start;
            let mut read_idx: usize = 0;

            cigar_per_base
                .zip(parse_md(md_tag))
                .for_each(|((r, q), md_base)| {
                    if r == 1 && q == 1 && md_base == Some(b'A') {
                        let read_base = sequence[read_idx];
                        let entry = map.entry((chr_id, ref_pos)).or_insert((0, 0));
                        match read_base {
                            b'G' => entry.0 += 1,
                            _ => entry.1 += 1,
                        }
                    }
                    ref_pos += r as u32;
                    read_idx += q as usize;
                });
        });

    let mut entries: Vec<((u8, u32), (u32, u32))> = map.into_iter().collect();
    entries.sort_by_key(|(key, _)| *key);

    let outfile = File::create(&cli.output_path)?;
    let mut writer = BufWriter::new(outfile);

    for ((chr_id, pos), (g, other)) in entries {
        writer.write_all(&[chr_id])?;
        writer.write_all(&pos.to_le_bytes())?;
        writer.write_all(&g.to_le_bytes())?;
        writer.write_all(&other.to_le_bytes())?;
    }
    writer.flush()?;

    Ok(())
}

fn extract_md_tag<'a>(fields: &'a [&str]) -> Option<&'a str> {
    fields
        .iter()
        .find(|field| field.starts_with("MD:Z:"))
        .map(|field| &field[5..])
}

fn parse_cigar(cigar: &str) -> impl Iterator<Item = (u32, u8, u8)> {
    cigar
        .split_inclusive(|c: char| c.is_ascii_alphabetic() || c == '=')
        .map(|pattern| {
            let split_pos = pattern.len() - 1;
            let (len, op) = pattern.split_at(split_pos);

            let len = len.parse::<u32>().expect("invalid CIGAR length");
            let op = op.as_bytes()[0] as usize;

            let (r, q) = CIGAR_TABLE[op];
            (len, r, q)
        })
}

enum MdToken {
    Matches(std::iter::Take<std::iter::Repeat<Option<u8>>>),
    Deletion,
    Mismatch(std::iter::Once<Option<u8>>),
}

impl Iterator for MdToken {
    type Item = Option<u8>;
    fn next(&mut self) -> Option<Option<u8>> {
        match self {
            MdToken::Matches(it) => it.next(),
            MdToken::Deletion => None,
            MdToken::Mismatch(it) => it.next(),
        }
    }
}

fn parse_md(md: &str) -> impl Iterator<Item = Option<u8>> {
    md.as_bytes()
        .chunk_by(|a, b| a.is_ascii_digit() == b.is_ascii_digit())
        .flat_map(|chunk| {
            if chunk[0].is_ascii_digit() {
                let n: usize = std::str::from_utf8(chunk).unwrap().parse().unwrap();
                MdToken::Matches(std::iter::repeat(None).take(n))
            } else if chunk[0] == b'^' {
                MdToken::Deletion
            } else {
                MdToken::Mismatch(std::iter::once(Some(chunk[0])))
            }
        })
}
