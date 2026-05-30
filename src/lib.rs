//! Extend BED intervals by N bp on each side — bedtools slop equivalent.
//!
//! Coordinates are clamped to [0, chrom_size]; unknown chromosomes pass through
//! unchanged with a stderr warning (matches bedtools).

use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub fn read_genome(path: &Path) -> Result<HashMap<String, u64>> {
    let data = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut map = HashMap::new();
    for (lineno, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.splitn(2, '\t');
        let chrom = fields.next().unwrap_or("").to_owned();
        let size_str = fields.next().unwrap_or("");
        let size: u64 = size_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!(
                "genome file {}: line {}: bad size {:?}",
                path.display(),
                lineno + 1,
                size_str
            ))
        })?;
        map.insert(chrom, size);
    }
    Ok(map)
}

pub struct SlopConfig {
    pub left: u64,
    pub right: u64,
    /// When true, left/right are percentages of interval length rather than bp.
    pub pct: bool,
}

impl SlopConfig {
    pub fn symmetric(b: u64) -> Self {
        Self {
            left: b,
            right: b,
            pct: false,
        }
    }
}

pub fn slop(
    input: &Path,
    genome: &HashMap<String, u64>,
    cfg: &SlopConfig,
    output: &mut dyn Write,
) -> Result<()> {
    let file = std::fs::File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    slop_reader(BufReader::new(file), genome, cfg, output)
}

pub fn slop_stdin(
    genome: &HashMap<String, u64>,
    cfg: &SlopConfig,
    output: &mut dyn Write,
) -> Result<()> {
    slop_reader(BufReader::new(io::stdin()), genome, cfg, output)
}

fn slop_reader<R: io::Read>(
    reader: BufReader<R>,
    genome: &HashMap<String, u64>,
    cfg: &SlopConfig,
    output: &mut dyn Write,
) -> Result<()> {
    let mut out = BufWriter::new(output);
    let mut stderr = io::stderr().lock();

    for (lineno_0, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let bytes = line.as_bytes();

        if bytes.is_empty()
            || bytes[0] == b'#'
            || bytes.starts_with(b"track")
            || bytes.starts_with(b"browser")
        {
            out.write_all(bytes).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }

        let lineno = lineno_0 + 1;
        let mut fields = line.splitn(4, '\t');
        let chrom = fields.next().unwrap_or("");
        let start_str = fields.next().unwrap_or("");
        let end_str = fields.next().unwrap_or("");
        let rest = fields.next().unwrap_or("");

        let start: u64 = start_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad start {start_str:?}"))
        })?;
        let end: u64 = end_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad end {end_str:?}"))
        })?;

        let (left_bp, right_bp) = if cfg.pct {
            let len = end.saturating_sub(start) as f64;
            (
                (cfg.left as f64 / 100.0 * len).round() as u64,
                (cfg.right as f64 / 100.0 * len).round() as u64,
            )
        } else {
            (cfg.left, cfg.right)
        };

        let chrom_size = match genome.get(chrom) {
            Some(&s) => s,
            None => {
                writeln!(stderr, "Warning: chromosome {chrom:?} not in genome file — interval passed through unchanged").ok();
                // Pass through unchanged.
                out.write_all(bytes).map_err(RsomicsError::Io)?;
                out.write_all(b"\n").map_err(RsomicsError::Io)?;
                continue;
            }
        };

        let new_start = start.saturating_sub(left_bp);
        let new_end = (end + right_bp).min(chrom_size);

        write!(out, "{chrom}\t{new_start}\t{new_end}").map_err(RsomicsError::Io)?;
        if !rest.is_empty() {
            write!(out, "\t{rest}").map_err(RsomicsError::Io)?;
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
