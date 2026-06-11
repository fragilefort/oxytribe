use std::fs::File;
use std::io::{BufRead, BufReader};

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
