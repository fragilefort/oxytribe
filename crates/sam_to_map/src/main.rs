use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let f = File::open("sample.sam").expect("Couldn't open sam file");
    let reader = BufReader::new(f);
    let mut lines_iter = reader.lines().map(|l| l.unwrap());
    lines_iter.next().filter(|line| line.starts_with("@"));
    Ok(())
}
