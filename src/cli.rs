use clap::Parser;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// a FASTA-formatted file
    #[arg(value_name = "FILE", required = true)]
    fasta: String,

    /// a list of regions to extract in SAMtools region format (chr1:1-1000, chr1);
    /// a negative sign in front of a region causes the extracted region to be reverse complemented
    #[arg(value_name = "FILE", required = true)]
    regions: String,

    /// output to this location (default is stdout)
    #[arg(short, long, value_name = "FILE")]
    output: Option<String>,

    /// default output is individual regions/contigs; this options outputs a single merged contig
    #[arg(short, long)]
    merge_contigs: bool,

    /// name of the single merged contig (default is regions filename without extension)
    #[arg(short, requires = "merge_contigs")]
    contig_name: Option<String>,

    /// insert gaps of this length between sequences
    #[arg(short, requires = "merge_contigs", default_value_t = 0)]
    gap_size: usize,
}

impl Cli {
    pub fn get_input(&self) -> (String, String) {
        (self.fasta.clone(), self.regions.clone())
    }

    pub fn get_output(&self) -> (Option<String>, bool, Option<String>, usize) {
        (
            self.output.clone(),
            self.merge_contigs,
            self.contig_name.clone(),
            self.gap_size,
        )
    }
}
