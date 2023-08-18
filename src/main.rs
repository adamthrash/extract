use std::collections::HashMap;
use std::fs::{read_to_string, File};
use std::path::Path;

use clap::Parser;
use noodles::core::Position;
use noodles::fasta::{self as fasta, fai, Record};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// a FASTA-formatted file
    #[arg(value_name = "FILE")]
    fasta: String,
    /// a list of regions to extract in SAMtools region format (chr1:1-1000, chr1);
    /// a negative sign in front of a region causes the extracted region to be reversed
    #[arg(value_name = "FILE")]
    regions: String,
    /// by default, regions (such as chr1:1-10 and chr1:1-20) are not merged;
    /// this option merges them into a single sequence (chr1)
    #[arg(short, long)]
    merge_regions: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    let fasta_file = args.fasta;
    let region_file = args.regions;

    let mut reader = if std::path::Path::new(&format! {"{fasta_file}.fai"}).exists() {
        fasta::indexed_reader::Builder::default().build_from_path(fasta_file)?
    } else {
        let file = File::create(format! {"{fasta_file}.fai"})?;
        let index = fasta::index(fasta_file.clone())?;
        let mut writer = fai::Writer::new(file);
        writer.write_index(&index)?;
        fasta::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(fasta_file)?
    };

    let mut ordered_sequences: Vec<String> = Vec::new();
    let mut sequences: HashMap<String, Record> = HashMap::new();

    for raw_region in read_to_string(region_file.clone()).unwrap().lines() {
        let mut reverse = false;
        let mut raw_region = raw_region.to_string();
        if raw_region.starts_with('-') {
            reverse = true;
            raw_region = raw_region[1..].to_string();
        }

        let region = raw_region.parse()?;
        let mut record = reader.query(&region)?;

        if reverse {
            let definition = fasta::record::Definition::new(record.name(), None);
            let start = Position::try_from(1)?;
            let end = Position::try_from(record.sequence().len())?;
            let mut sequence = record
                .sequence()
                .get(start..=end)
                .expect("could not get sequence")
                .to_vec();
            sequence.reverse();
            record = fasta::Record::new(definition, sequence.into());
        }

        let name = if args.merge_regions {
            if let Some(split) = record.name().to_string().split_once(':') {
                split.0.to_string()
            } else {
                record.name().to_string()
            }
        } else {
            record.name().to_string()
        };

        if !ordered_sequences.contains(&name) {
            ordered_sequences.push(name.clone());
        }
        let definition = fasta::record::Definition::new(name.clone(), None);
        let start = Position::try_from(1)?;
        let end = Position::try_from(record.sequence().len())?;
        let sequence = record
            .sequence()
            .get(start..=end)
            .expect("could not get sequence")
            .to_vec();
        record = fasta::Record::new(definition, sequence.into());

        sequences
            .entry(name)
            .and_modify(|existing_record| {
                let start = Position::try_from(1).expect("could not get position");
                let end = Position::try_from(existing_record.sequence().len())
                    .expect("could not get position");
                let mut sequence = existing_record
                    .sequence()
                    .get(start..=end)
                    .expect("could not get sequence")
                    .to_vec();
                let start = Position::try_from(1).expect("could not get position");
                let end =
                    Position::try_from(record.sequence().len()).expect("could not get position");
                let mut extended_sequence = record
                    .sequence()
                    .get(start..=end)
                    .expect("could not get sequence")
                    .to_vec();
                sequence.append(&mut extended_sequence);
                let record = fasta::Record::new(
                    fasta::record::Definition::new(existing_record.name(), None),
                    sequence.into(),
                );
                *existing_record = record;
            })
            .or_insert(record.clone());
    }

    let definition = fasta::record::Definition::new(
        Path::new(&region_file)
            .file_stem()
            .unwrap()
            .to_str()
            .expect("could not get str"),
        None,
    );
    let mut scaffold = fasta::Record::new(definition, Vec::new().into());

    for key in ordered_sequences {
        let start = Position::try_from(1).expect("could not get position");
        let mut sequence = if let Ok(end) = Position::try_from(scaffold.sequence().len()) {
            scaffold
                .sequence()
                .get(start..=end)
                .expect("could not get sequence")
                .to_vec()
        } else {
            Vec::new()
        };
        let record = &sequences.get(&key).expect("could not get key");
        let start = Position::try_from(1).expect("could not get position");
        let end = Position::try_from(record.sequence().len()).expect("could not get position");
        let mut extended_sequence = record
            .sequence()
            .get(start..=end)
            .expect("could not get sequence")
            .to_vec();
        sequence.append(&mut extended_sequence);
        scaffold = fasta::Record::new(
            fasta::record::Definition::new(scaffold.name(), None),
            sequence.into(),
        );
    }

    let file = File::create(format!(
        "{}.fasta",
        Path::new(&region_file)
            .file_stem()
            .unwrap()
            .to_str()
            .expect("could not get str")
    ))?;
    let mut writer = fasta::Writer::new(file);
    writer.write_record(&scaffold)?;

    Ok(())
}
