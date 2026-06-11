use std::fs::File;
use std::io::BufReader;

fn main() -> std::io::Result<()> {
    let f = File::open("sample.sam").expect("Couldn't open sam file");
    let mut reader = BufReader::new(f);

    Ok(())
}
