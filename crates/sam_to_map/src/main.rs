use either::Either;
use std::collections::HashMap;
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
    let mut map: HashMap<u32, (u32, u32)> = HashMap::with_capacity(1_000_000);
    let f = File::open("sample.sam").expect("Couldn't open sam file");
    let reader = BufReader::new(f);
    reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|line| !line.starts_with("@"))
        .for_each(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let start = fields[3].parse::<u32>().unwrap();
            let cigar = fields[5];
            let sequence = fields[9].as_bytes();
            let md_tag = extract_md_tag(&fields).unwrap_or("");

            parse_cigar(cigar).zip(parse_md(md_tag))
            // emit a2g and insert into map
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

fn parse_md(md: &str) -> impl Iterator<Item = Option<u8>> {
    md.as_bytes()
        .chunk_by(|a, b| a.is_ascii_digit() == b.is_ascii_digit())
        .flat_map(|chunk| {
            if chunk[0].is_ascii_digit() {
                let n: usize = std::str::from_utf8(chunk).unwrap().parse().unwrap();
                Either::Left(std::iter::repeat(None).take(n))
            } else {
                Either::Right(std::iter::once(Some(chunk[0])))
            }
        })
}
