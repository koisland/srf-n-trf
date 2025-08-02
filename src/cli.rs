use std::{f32, path::PathBuf};

use clap::Parser;

/// Script to take `srf` and `trf` output and produce a bed file with only regions corresponding monomers of a given periodicity.
#[derive(Debug, Parser)]
#[command(version, about, long_about = None)]
pub struct Cli {
    #[arg(short, long)]
    /// PAF file of alignment of assembly as query and `srf` enlonged motifs as target.
    /// Requires `cg` extended cigar string. With `minimap2`, use `--eqx`.
    pub paf: PathBuf,
    /// `trf` monomers TSV file on `srf` monomers with columns:
    /// `chrom (query), motif (target), st, end, period, copyNum, fracMatch, fracGap, score, entropy, pattern`
    #[arg(short, long)]
    pub monomers: PathBuf,
    /// Output BED9 file with columns:
    /// `chrom, st, end, comma-delimited_monomers, 0, strand, st, end, '0,0,0'`
    #[arg(short, long)]
    pub outfile: Option<PathBuf>,
    /// Monomer size in base pairs to search for.
    #[arg(short, long, default_values_t = [170, 340, 42])]
    pub sizes: Vec<u32>,
    /// Percent difference in monomer period allowed.
    /// ex. `0.02` results in valid periods for `170`: `167 < 170 < 173`
    #[arg(short, long, default_value_t = 0.02)]
    pub diff: f32,
}
