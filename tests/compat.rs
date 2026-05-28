use std::path::Path;
use std::process::Command;

use rsomics_bed_slop::{SlopConfig, read_genome, slop};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn basic_slop_both_sides() {
    let input = golden("input.bed");
    let genome_path = golden("genome.txt");
    let genome = read_genome(&genome_path).unwrap();
    let cfg = SlopConfig::symmetric(50);
    let mut out = Vec::new();
    slop(&input, &genome, &cfg, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    // chr1:100-200 + 50 each side → chr1:50-250
    assert!(result.contains("chr1\t50\t250"), "A not extended: {result}");
    // chr1:900-950 + 50 each side → chr1:850-1000 (clamped at chrom size)
    assert!(
        result.contains("chr1\t850\t1000"),
        "B not clamped: {result}"
    );
    // chr2:50-150 + 50 each side → chr2:0-200
    assert!(result.contains("chr2\t0\t200"), "C not extended: {result}");
}

#[test]
fn start_clamps_at_zero() {
    let input = golden("input.bed");
    let genome_path = golden("genome.txt");
    let genome = read_genome(&genome_path).unwrap();
    let cfg = SlopConfig::symmetric(200);
    let mut out = Vec::new();
    slop(&input, &genome, &cfg, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    // chr1:100-200 - 200 = 0 (clamped)
    let lines: Vec<&str> = result.lines().filter(|l| l.starts_with("chr1\t")).collect();
    let starts: Vec<u64> = lines
        .iter()
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect();
    for s in &starts {
        assert!(*s == 0, "start should be clamped to 0, got: {s}");
    }
}

#[test]
fn extra_columns_preserved() {
    let input = golden("input.bed");
    let genome_path = golden("genome.txt");
    let genome = read_genome(&genome_path).unwrap();
    let cfg = SlopConfig::symmetric(10);
    let mut out = Vec::new();
    slop(&input, &genome, &cfg, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    // BED6 name/score/strand columns should survive
    assert!(result.contains("regionA"), "name column lost: {result}");
    assert!(result.contains("regionB"), "name column lost: {result}");
}

#[test]
fn bedtools_compat() {
    let bedtools = Command::new("bedtools").arg("--version").output();
    if bedtools.is_err() || !bedtools.unwrap().status.success() {
        eprintln!("bedtools not available — skipping compat test");
        return;
    }

    let input = golden("input.bed");
    let genome_path = golden("genome.txt");
    let genome = read_genome(&genome_path).unwrap();
    let cfg = SlopConfig::symmetric(50);

    let mut ours = Vec::new();
    slop(&input, &genome, &cfg, &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    let bt = Command::new("bedtools")
        .args(["slop", "-i"])
        .arg(&input)
        .arg("-g")
        .arg(&genome_path)
        .args(["-b", "50"])
        .output()
        .expect("bedtools slop failed");
    let bt_str = String::from_utf8(bt.stdout).unwrap();

    let mut ours_lines: Vec<&str> = ours_str.lines().filter(|l| !l.is_empty()).collect();
    let mut bt_lines: Vec<&str> = bt_str.lines().filter(|l| !l.is_empty()).collect();
    ours_lines.sort_unstable();
    bt_lines.sort_unstable();

    assert_eq!(ours_lines, bt_lines, "output differs from bedtools slop");
}
