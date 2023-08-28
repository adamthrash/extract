use clap::Parser;

mod cli;
mod sequences;

use anyhow::Result;
use cli::Cli;
use sequences::Sequences;

fn main() -> Result<()> {
    // Parse CLI arguments
    let args = Cli::parse();
    let (fasta_file, region_file) = args.get_input();
    let (output_location, merge, contig_name, gap_size) = args.get_output();

    // Create Sequences struct; extract sequences; write output.
    let mut sequences = Sequences::new(&fasta_file, &region_file)?;
    sequences.extract()?;
    sequences.write(output_location, merge, contig_name, gap_size)?;
    Ok(())
}
