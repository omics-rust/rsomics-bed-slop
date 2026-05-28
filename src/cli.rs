use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_slop::{SlopConfig, read_genome, slop, slop_stdin};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-slop", disable_help_flag = true)]
pub struct Cli {
    /// Input BED (default: stdin)
    input: Option<PathBuf>,
    /// Chromosome sizes file (required; two-column chrom\tsize TSV)
    #[arg(short = 'g', long, value_name = "FILE")]
    genome: PathBuf,
    /// Output BED (default: stdout)
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    /// Extend both sides by N bp (shorthand for -l N -r N)
    #[arg(short = 'b', long, conflicts_with_all = ["left", "right"])]
    both: Option<u64>,
    /// Extend left (5′) side by N bp
    #[arg(short = 'l', long)]
    left: Option<u64>,
    /// Extend right (3′) side by N bp
    #[arg(short = 'r', long)]
    right: Option<u64>,
    /// Treat -l/-r/-b as percentage of interval length instead of bp
    #[arg(long, default_value = "false")]
    pct: bool,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let cfg = if let Some(b) = self.both {
            SlopConfig {
                left: b,
                right: b,
                pct: self.pct,
            }
        } else {
            let left = self.left.unwrap_or(0);
            let right = self.right.unwrap_or(0);
            if self.left.is_none() && self.right.is_none() {
                return Err(RsomicsError::InvalidInput(
                    "one of -b/--both, -l/--left, or -r/--right is required".into(),
                ));
            }
            SlopConfig {
                left,
                right,
                pct: self.pct,
            }
        };

        let genome = read_genome(&self.genome)?;

        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };

        match self.input {
            Some(ref p) => slop(p.as_path(), &genome, &cfg, out),
            None => slop_stdin(&genome, &cfg, out),
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Extend BED intervals by N bp on each side, clamped to chrom bounds (bedtools slop).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -g <genome> [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<FILE>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Chromosome sizes file (chrom\\tsize TSV)",
                why_default: None,
            },
            FlagSpec {
                short: Some('b'),
                long: "both",
                aliases: &[],
                value: Some("<INT>"),
                type_hint: Some("u64"),
                required: false,
                default: None,
                description: "Extend both sides by N bp",
                why_default: None,
            },
            FlagSpec {
                short: Some('l'),
                long: "left",
                aliases: &[],
                value: Some("<INT>"),
                type_hint: Some("u64"),
                required: false,
                default: None,
                description: "Extend left (5′) side by N bp",
                why_default: None,
            },
            FlagSpec {
                short: Some('r'),
                long: "right",
                aliases: &[],
                value: Some("<INT>"),
                type_hint: Some("u64"),
                required: false,
                default: None,
                description: "Extend right (3′) side by N bp",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "pct",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Treat extension as percent of interval length",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output BED path",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Extend each interval by 100 bp on both sides",
            command: "rsomics-bed-slop -i input.bed -g chrom.sizes -b 100",
        },
        Example {
            description: "Extend 50 bp upstream, 200 bp downstream",
            command: "rsomics-bed-slop -i input.bed -g chrom.sizes -l 50 -r 200",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
