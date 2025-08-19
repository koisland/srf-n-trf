use std::{f32, path::PathBuf};

use clap::{Parser, Subcommand};

/// Script to take `srf` and `trf` output and produce a bed file with only regions corresponding monomers of a given periodicity.
#[derive(Debug, Parser)]
#[command(version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    Monomers {
        #[arg(short, long)]
        /// PAF file of alignment of assembly as query and `srf` enlonged motifs as target.
        /// Requires `cg` extended cigar string. With `minimap2`, use `--eqx`.
        paf: PathBuf,
        /// `trf` monomers TSV file on `srf` monomers with columns:
        /// `chrom (query), motif (target), st, end, period, copyNum, fracMatch, fracGap, score, entropy, pattern`
        #[arg(short, long)]
        monomers: PathBuf,
        /// Output BED9 file with columns:
        /// `chrom, st, end, comma-delimited_monomers, 0, strand, st, end, '0,0,0'`
        #[arg(short, long)]
        outfile: Option<PathBuf>,
        /// Monomer size in base pairs to search for.
        #[arg(short, long, default_values_t = [170, 340, 510, 680, 850, 1020, 42], num_args = 1..)]
        sizes: Vec<u32>,
        /// Percent difference in monomer period length allowed.
        /// ex. `0.02` results in valid periods for `170`: `167 < 170 < 173`
        #[arg(short, long, default_value_t = 0.02)]
        diff: f32,
        /// Maximum gap-compressed sequence divergence between aligned motif and region.
        #[arg(short, long, default_value_t = 0.2)]
        max_seq_div: f64,
    },
    Motifs {
        #[arg(short, long)]
        /// Fasta file of srf detected motifs.
        fa: PathBuf,
        /// `trf` monomers TSV file on `srf` monomers with columns:
        /// `chrom (query), motif (target), st, end, period, copyNum, fracMatch, fracGap, score, entropy, pattern`
        #[arg(short, long)]
        monomers: PathBuf,
        /// Output fasta file filtered to only motifs composed of monomers of given size.
        #[arg(short, long)]
        outfile: Option<PathBuf>,
        /// Monomer size in base pairs to search for.
        #[arg(short, long, default_values_t = [170, 340, 510, 680, 850, 1020], num_args = 1..)]
        sizes: Vec<u32>,
        /// Percent difference in monomer period length allowed.
        /// ex. `0.02` results in valid periods for `170`: `167 < 170 < 173`
        #[arg(short, long, default_value_t = 0.02)]
        diff: f32,
        /// Require all monomers to be within size range.
        #[arg(long, action)]
        require_all: bool,
    },
    Regions {
        /// Bed file from extract command.
        #[arg(short, long)]
        bed: PathBuf,
        /// Output BED9 file with columns:
        /// `chrom, st, end, comma-delimited_monomers, 0, strand, st, end, '0,0,0'`
        #[arg(short, long)]
        outfile: Option<PathBuf>,
        /// Distance to merge in base pairs.
        #[arg(short, long, default_value_t = 100_000)]
        dst: u32,
        /// Minimum length in base pairs.
        #[arg(short, long, default_value_t = 30_000)]
        min_len: u32,
        /// Required monomers in merged blocks. Merges iff one of these monomer periods is in block.
        /// Also filters out monomers not within this period.
        #[arg(short, long, default_values_t = [170, 340, 510, 680, 850, 1020], num_args = 1..)]
        sizes: Vec<u32>,
        /// Difference in required monomer size.
        #[arg(long, default_value_t = 0.02)]
        diff: f32,
    },
}
