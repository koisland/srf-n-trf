use std::{
    collections::{HashSet, VecDeque},
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write, stdin, stdout},
};

use clap::Parser;
use eyre::{ContextCompat, bail};
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

pub fn create_monomer_range(sizes: &[u32], diff: f32) -> Lapper<u32, ()> {
    Lapper::new(
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
    )
}

/// `writeln!()` but handles broken pipes.
/// * https://stackoverflow.com/a/65760807
macro_rules! writeln_w_bp {
    ($dst:expr, $($arg:tt)*) => {
        if let Err(err) = writeln!($dst, $($arg)*) {
            if err.kind() != std::io::ErrorKind::BrokenPipe {
                let _ = writeln!(std::io::stderr(), "{err:?}");
            }
        }
    };
}

fn main() -> eyre::Result<()> {
    let cli = Cli::parse();
    eprintln!("Running command:\n{:#?}", &cli.command);

    match cli.command {
        Command::Monomers {
            paf,
            monomers,
            outfile,
            sizes,
            diff,
            max_seq_div,
        } => {
            let monomers = read_trf_monomers(monomers)?;
            let reader = Reader::from_path(paf)?;
            // Inteval tree of allowed period ranges.
            let monomer_period_range: Lapper<u32, ()> = create_monomer_range(&sizes, diff);
            let monomer_periods: HashSet<u32> = HashSet::from_iter(sizes);
            let mut writer = if let Some(outfile) = outfile {
                Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
            } else {
                Box::new(BufWriter::new(stdout().lock())) as Box<dyn Write>
            };

            let Some(min_monomer_period) = monomer_periods.iter().min().cloned() else {
                bail!("No monomer periods provided.")
            };

            eprintln!(
                "Using monomer periodicity range:\n{:#?}",
                monomer_period_range.intervals
            );
            let null_lapper = Lapper::new(vec![]);

            for rec in reader
                .into_records()
                .flatten()
                .sorted_by(|a, b| a.query_start().cmp(&b.query_start()))
            {
                let target_tr_chrom_monomers = monomers
                    .get(rec.query_name())
                    .and_then(|mp| mp.get(rec.target_name()))
                    .unwrap_or(&null_lapper);
                let target_len = rec.target_len() as i32;
                let aln_len = rec.alignment_block_len() as i32;
                let aln_itv_diff = target_len.abs_diff(aln_len);
                let aln_rpt_len_perc_diff = aln_itv_diff as f32 / rec.target_len() as f32;

                // If rec is within x% difference in length. Use gap-comprssed identity rather than overlap to find divergent and monomeric HORs.
                // Will not return individual monomer positions but entire region.
                if aln_rpt_len_perc_diff < diff
                    && rec.de().map(|de| *de < max_seq_div).unwrap_or_default()
                {
                    let mut monomers = target_tr_chrom_monomers
                        .iter()
                        .filter_map(|m| {
                            (monomer_period_range.count(m.val.trf_period, m.val.trf_period) > 0)
                                .then_some(&m.val.trf_monomer)
                        })
                        .join(",");

                    // Allow if motif found is within range even if doesn't haven any monomers.
                    if monomer_period_range
                        .count(rec.alignment_block_len(), rec.alignment_block_len())
                        > 0
                    {
                        monomers.push('.');
                    } else if monomers.is_empty() {
                        continue;
                    }
                    writeln_w_bp!(
                        &mut writer,
                        "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t0,0,0",
                        rec.query_name(),
                        rec.query_start(),
                        rec.query_end(),
                        monomers,
                        rec.strand(),
                        rec.query_start(),
                        rec.query_end(),
                    );
                    continue;
                }

                // Otherwise, search cigar string elements for monomers.
                let paired_itvs = get_aligned_paired_itvs(&rec, min_monomer_period)?;
                for (q_itv, t_itv) in paired_itvs {
                    let ovl = target_tr_chrom_monomers
                        .find(t_itv.start, t_itv.stop)
                        .collect_vec();

                    if ovl.is_empty() {
                        continue;
                    }
                    let q_itv_len = q_itv.stop - q_itv.start;

                    let monomers = ovl
                        .iter()
                        .filter_map(|o| {
                            let is_period_ovl =
                                monomer_period_range.count(o.val.trf_period, o.val.trf_period) > 0;
                            let at_least_period_size = o.val.trf_period <= q_itv_len;

                            (is_period_ovl && at_least_period_size).then_some(&o.val.trf_monomer)
                        })
                        .join(",");

                    if monomers.is_empty() {
                        continue;
                    }
                    writeln_w_bp!(
                        &mut writer,
                        "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t0,0,0",
                        rec.query_name(),
                        q_itv.start,
                        q_itv.stop,
                        monomers,
                        rec.strand(),
                        q_itv.start,
                        q_itv.stop,
                    );
                }
            }
        }
        Command::Motifs {
            fa,
            monomers,
            outfile,
            sizes,
            diff,
            require_all,
        } => {
            let reader = if fa != OsStr::new("-") {
                Box::new(BufReader::new(File::open(fa)?)) as Box<dyn BufRead>
            } else {
                Box::new(BufReader::new(stdin().lock()))
            };
            let monomers = read_trf_monomers(monomers)?;
            let (_, monomers) = monomers.into_iter().next().context("No monomers.")?;

            let mut writer = if let Some(outfile) = outfile {
                Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
            } else {
                Box::new(BufWriter::new(stdout().lock())) as Box<dyn Write>
            };
            let monomer_period_range: Lapper<u32, ()> = create_monomer_range(&sizes, diff);

            let mut curr_rec: Option<String> = None;
            for line in reader.lines().map_while(Result::ok) {
                if line.starts_with('>') {
                    curr_rec = line.trim().strip_prefix('>').and_then(|name| {
                        name.split_once(' ')
                            .map(|(name, _comments)| name.to_owned())
                    });
                    continue;
                }
                if let Some((rec_name, rec_monomers)) = curr_rec
                    .as_ref()
                    .and_then(|rec_name| monomers.get(rec_name).map(|mons| (rec_name, mons)))
                {
                    let is_valid_motif = if require_all {
                        rec_monomers.iter().all(|mon| {
                            monomer_period_range.count(mon.val.trf_period, mon.val.trf_period) > 0
                        })
                    } else {
                        rec_monomers.iter().any(|mon| {
                            monomer_period_range.count(mon.val.trf_period, mon.val.trf_period) > 0
                        })
                    };
                    if is_valid_motif {
                        writeln_w_bp!(&mut writer, ">{rec_name}");
                        writeln_w_bp!(&mut writer, "{line}");
                    }
                }
            }
        }
        Command::Regions {
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
                Box::new(BufReader::new(stdin().lock()))
            };
            let mut writer = if let Some(outfile) = outfile {
                Box::new(BufWriter::new(File::create(outfile)?)) as Box<dyn Write>
            } else {
                Box::new(BufWriter::new(stdout().lock())) as Box<dyn Write>
            };
            let monomer_period_range: Lapper<u32, ()> = create_monomer_range(&sizes, diff);
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

            // Merge intervals.
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
                // Must be same name and within distance.
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
                    // Filter monomers that don't fall within monomer period range.
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
                writeln_w_bp!(
                    &mut writer,
                    "{chrom}\t{st}\t{end}\t{monomers}\t0\t.\t{st}\t{end}\t0,0,0"
                );
            }
        }
    }

    Ok(())
}
