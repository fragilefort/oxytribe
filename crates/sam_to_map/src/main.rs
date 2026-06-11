use std::fs::File;
use std::io::{BufRead, BufReader};

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

fn main() -> std::io::Result<()> {
    let f = File::open("sample.sam").expect("Couldn't open sam file");
    let reader = BufReader::new(f);
    let lines_iter = reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|line| !line.starts_with("@"))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let start = fields[3].parse::<u32>();
            let cigar = fields[5];
            let sequence = fields[9];
            let md_tag = extract_md_tag(&fields);
        });

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
        .split_inclusive(|c: char| c.is_ascii_alphabetic())
        .map(|pattern| {
            let split_pos = pattern.len() - 1;
            let (len, op) = pattern.split_at(split_pos);

            let len = len.parse::<u32>().expect("invalid CIGAR length");
            let op = op.chars().next().expect("missing CIGAR op") as usize;

            let (r, q) = CIGAR_TABLE[op];
            (len, r, q)
        })
}
