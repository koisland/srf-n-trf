use std::{
    collections::HashSet,
    fs::File,
    io::{BufWriter, Write, stdout},
};

use clap::Parser;
use eyre::bail;
use itertools::Itertools;
use paf::Reader;
use rust_lapper::{Interval, Lapper};

mod cigar;
mod cli;
mod io;

use crate::{cigar::get_aligned_paired_itvs, cli::Cli, io::read_trf_monomers};

fn main() -> eyre::Result<()> {
    let cli = Cli::parse();
    eprintln!("Using config:\n{:#?}", &cli);

    let monomers = read_trf_monomers(cli.monomers)?;
    let reader = Reader::from_path(cli.paf)?;
    let monomer_periods: HashSet<u32> = HashSet::from_iter(cli.sizes);
    let perc_diff = cli.diff;
    let mut writer = if let Some(outfile) = cli.outfile {
        Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
    } else {
        Box::new(BufWriter::new(stdout())) as Box<dyn Write>
    };

    let Some(min_monomer_period) = monomer_periods.iter().min().cloned() else {
        bail!("No monomer periods provided.")
    };

    // Inteval tree of allowed period ranges.
    let monomer_period_range: Lapper<u32, ()> = Lapper::new(
        monomer_periods
            .iter()
            .map(|period| {
                let allowed_diff = *period as f32 * perc_diff;
                let lower_bound = *period as f32 - allowed_diff;
                let upper_bound = *period as f32 + allowed_diff;
                Interval {
                    start: lower_bound as u32,
                    stop: upper_bound as u32 + 1,
                    val: (),
                }
            })
            .collect(),
    );
    eprintln!(
        "Using monomer periodicity range:\n{:#?}",
        monomer_period_range.intervals
    );

    for rec in reader
        .into_records()
        .flatten()
        .filter(|rec| rec.tp().eq(&Some(&'P')))
        .sorted_by(|a, b| a.query_start().cmp(&b.query_start()))
    {
        let paired_itvs = get_aligned_paired_itvs(&rec, min_monomer_period)?;

        let Some(target_chrom_monomers) = monomers.get(rec.query_name()) else {
            continue;
        };
        for (q_itv, t_itv) in paired_itvs {
            let Some(target_tr_chrom_monomers) = target_chrom_monomers.get(rec.target_name())
            else {
                continue;
            };
            let ovl = target_tr_chrom_monomers
                .find(t_itv.start, t_itv.stop)
                .collect_vec();

            if ovl.is_empty() {
                continue;
            }

            let monomers = ovl
                .iter()
                .filter_map(|o| {
                    (monomer_period_range.count(o.val.trf_period, o.val.trf_period) > 0)
                        .then_some(&o.val.trf_monomer)
                })
                .join(",");
            if monomers.is_empty() {
                continue;
            }
            writeln!(
                &mut writer,
                "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t0,0,0",
                rec.query_name(),
                q_itv.start,
                q_itv.stop,
                monomers,
                rec.strand(),
                q_itv.start,
                q_itv.stop,
            )?;
        }
    }

    Ok(())
}
