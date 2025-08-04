use std::{
    collections::{HashSet, VecDeque},
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write, stdin, stdout},
};

use clap::Parser;
use eyre::bail;
use itertools::Itertools;
use paf::Reader;
use rust_lapper::{Interval, Lapper};

mod cigar;
mod cli;
mod io;

use crate::{
    cigar::get_aligned_paired_itvs,
    cli::{Cli, Command},
    io::read_trf_monomers,
};

fn main() -> eyre::Result<()> {
    let cli = Cli::parse();
    eprintln!("Running command:\n{:#?}", &cli.command);

    match cli.command {
        Command::Extract {
            paf,
            monomers,
            outfile,
            sizes,
            diff,
        } => {
            let monomers = read_trf_monomers(monomers)?;
            let reader = Reader::from_path(paf)?;
            let monomer_periods: HashSet<u32> = HashSet::from_iter(sizes);
            let mut writer = if let Some(outfile) = outfile {
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
                        let allowed_diff = *period as f32 * diff;
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
                    let Some(target_tr_chrom_monomers) =
                        target_chrom_monomers.get(rec.target_name())
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
        }
        Command::Merge {
            bed,
            outfile,
            dst,
            min_len,
            sizes,
            diff,
        } => {
            let reader = if bed != OsStr::new("-") {
                Box::new(BufReader::new(File::open(bed)?)) as Box<dyn BufRead>
            } else {
                Box::new(BufReader::new(stdin()))
            };
            let mut writer = if let Some(outfile) = outfile {
                Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
            } else {
                Box::new(BufWriter::new(stdout())) as Box<dyn Write>
            };
            let monomer_period_range: Lapper<u32, ()> = Lapper::new(
                sizes
                    .iter()
                    .map(|period| {
                        let allowed_diff = *period as f32 * diff;
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

            let mut final_intervals: Vec<(String, u32, u32, HashSet<String>)> = vec![];
            let mut intervals: VecDeque<(String, u32, u32, HashSet<String>)> = reader
                .lines()
                .map_while(Result::ok)
                .map(|line| {
                    let (chrom, st, end, monomers, _, _, _, _, _) =
                        line.split('\t').collect_tuple().unwrap();
                    (
                        chrom.to_owned(),
                        st.parse::<u32>().unwrap(),
                        end.parse::<u32>().unwrap(),
                        monomers.split(',').map(|m| m.to_owned()).collect(),
                    )
                })
                .sorted_by(|a, b| (&a.0, a.1).cmp(&(&b.0, b.1)))
                .collect();

            while let Some(mut itv_1) = intervals.pop_front() {
                let Some(itv_2) = intervals.pop_front() else {
                    // Remove anything that isn't in required monomer period range.
                    itv_1.3.retain(|m| {
                        monomer_period_range.count(m.len() as u32, m.len() as u32) != 0
                    });

                    final_intervals.push(itv_1);
                    break;
                };
                let dst_between = itv_2.1.saturating_sub(itv_1.2);
                if dst_between <= dst && itv_1.0 == itv_2.0 {
                    intervals.push_front((
                        itv_1.0,
                        itv_1.1,
                        itv_2.2,
                        itv_1
                            .3
                            .union(&itv_2.3)
                            .cloned()
                            .collect::<HashSet<String>>(),
                    ));
                } else {
                    let final_itv_len = itv_1.2 - itv_1.1;
                    itv_1.3.retain(|m| {
                        monomer_period_range.count(m.len() as u32, m.len() as u32) != 0
                    });

                    if final_itv_len > min_len && !itv_1.3.is_empty() {
                        final_intervals.push(itv_1);
                    }
                    intervals.push_front(itv_2);
                }
            }
            for (chrom, st, end, monomers) in final_intervals {
                let monomers = monomers.iter().join(",");
                writeln!(
                    &mut writer,
                    "{chrom}\t{st}\t{end}\t{monomers}\t0\t.\t{st}\t{end}\t0,0,0"
                )?;
            }
        }
    }

    Ok(())
}
