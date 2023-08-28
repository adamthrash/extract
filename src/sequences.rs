use std::{
    collections::HashMap,
    fs::{read_to_string, File},
    io::{self, Write},
    path::Path,
    str,
};

use anyhow::Result;
use noodles::{
    core::{Position, Region},
    fasta::{self as fasta, fai, io::BufReadSeek, record::Sequence, IndexedReader, Record},
};

// The Sequences struct contains
// - the order in which sequences should be printed
// - the regions as parsed
// - the FASTA file reader
// - a list of regions and whether the region is reverse complemented
// - the file stem of the regions file
pub struct Sequences {
    pub order: Vec<String>,
    pub data: HashMap<String, Record>,
    reader: IndexedReader<Box<dyn BufReadSeek>>,
    regions: Vec<(Region, bool)>,
    regions_filename: String,
}

impl Sequences {
    // Creating a Sequences struct initializes a blank Vec and HashMap for
    // the order and data respectively. It initializes the reader and
    // parses the regions file.
    pub fn new(fasta_file: &str, region_file: &str) -> Result<Self> {
        Ok(Self {
            order: Vec::new(),
            data: HashMap::new(),
            reader: Self::get_reader(fasta_file)?,
            regions: Self::get_regions(region_file)?,
            regions_filename: Path::new(&region_file)
                .file_stem()
                .unwrap()
                .to_str()
                .expect("could not get str")
                .to_string(),
        })
    }

    // Extracting the regions in a Sequence struct iterates of the regions
    // data and reverse complements the extracted record if necessary.
    // The order and record are stored.
    pub fn extract(&mut self) -> Result<()> {
        for (region, reversed) in &self.regions {
            let mut record = self.reader.query(region)?;
            if *reversed {
                let definition = fasta::record::Definition::new(record.name(), None);
                let sequence: Sequence = record
                    .sequence()
                    .complement()
                    .rev()
                    .collect::<Result<_, _>>()?;
                record = fasta::Record::new(definition, sequence);
            }
            let record_name = record.name().to_string();
            self.order.push(record_name.clone());
            self.data.insert(record_name, record);
        }
        Ok(())
    }

    // Writing output from a Sequences struct checks
    // - whether the output location is a file or stdout
    // - whether all contigs or a single merged contig should be written
    // - what the name of the single merged contig should be
    // - whether the single merged contig should have gaps of a specific size
    pub fn write(
        &self,
        output_location: Option<String>,
        merge: bool,
        contig_name: Option<String>,
        gap_size: usize,
    ) -> Result<()> {
        // Get a Writer to stdout or a file.
        let mut writer: fasta::Writer<Box<dyn Write>> = match output_location {
            Some(path) => fasta::Writer::new(Box::new(File::create(path)?)),
            None => fasta::Writer::new(Box::new(io::stdout().lock())),
        };

        if !merge {
            // If the user didn't request a merged contig, write each contig.
            for key in &self.order {
                let record = &self.data.get(key).expect("could not get key");
                writer.write_record(record)?;
            }
        } else {
            // Create a gap if the user specified a gap size.
            let gap = if gap_size > 0 {
                Some("N".repeat(gap_size))
            } else {
                None
            };

            // Iterate over the sequence data in order, extracting the sequence data from
            // the record and converting it to &str. Store the sequence data in a Vec, and
            // add the gap sequence if it exists. The resulting Vec<&str> is flattened, and
            // the Vec of sequence data (and optional gaps) is concatenated.
            let sequences: String = self
                .order
                .iter()
                .flat_map(|sequence| {
                    let record = &self.data.get(sequence).expect("could not get key");
                    let start = Position::try_from(1).expect("could not get position");
                    let end = Position::try_from(record.sequence().len())
                        .expect("could not get position");
                    let mut sequence_data = vec![str::from_utf8(
                        record
                            .sequence()
                            .get(start..=end)
                            .expect("could not get sequence"),
                    )
                    .expect("could not convert sequence to String")];
                    if let Some(gap) = &gap {
                        sequence_data.push(gap);
                    }
                    sequence_data
                })
                .collect::<Vec<&str>>()
                .join("");

            // Select the contig name from either user input or the regions file's name.
            let contig_name = if let Some(contig_name) = contig_name {
                contig_name
            } else {
                self.regions_filename.clone()
            };

            // Create and write the record.
            let definition = fasta::record::Definition::new(contig_name, None);
            let record = fasta::Record::new(definition, sequences.as_bytes().to_vec().into());
            writer.write_record(&record)?;
        }
        Ok(())
    }

    // Return an IndexedReader, creating an index if one does not exist.
    fn get_reader(fasta_file: &str) -> Result<IndexedReader<Box<dyn BufReadSeek>>> {
        Ok(
            if std::path::Path::new(&format! {"{fasta_file}.fai"}).exists() {
                fasta::indexed_reader::Builder::default().build_from_path(fasta_file)?
            } else {
                let file = File::create(format! {"{fasta_file}.fai"})?;
                let index = fasta::index(fasta_file)?;
                let mut writer = fai::Writer::new(file);
                writer.write_index(&index)?;
                fasta::indexed_reader::Builder::default()
                    .set_index(index)
                    .build_from_path(fasta_file)?
            },
        )
    }

    // Parse each non-blank line in the regions file, noting whether
    // it should be reverse complemented.
    fn get_regions(region_file: &str) -> Result<Vec<(Region, bool)>> {
        Ok(read_to_string(region_file)
            .unwrap()
            .lines()
            .filter_map(|region| {
                if region.is_empty() {
                    None
                } else {
                    let mut reverse = false;
                    let mut region = region.to_string();
                    if region.starts_with('-') {
                        reverse = true;
                        region = region[1..].to_string();
                    }

                    if let Ok(region) = region.parse() {
                        Some((region, reverse))
                    } else {
                        None
                    }
                }
            })
            .collect())
    }
}
