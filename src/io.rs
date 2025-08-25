use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use rust_lapper::{Interval, Lapper};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Monomer {
    pub srf_repeat: String,
    pub trf_monomer: String,
    pub trf_period: u32,
    pub trf_copy_num: OrderedFloat<f32>,
}

type MotifMonomers = HashMap<String, Lapper<u32, Monomer>>;

/*
INPUT_TRF_COLS = (
    "chrom",
    "motif",
    "st",
    "end",
    "period",
    "copyNum",
    "fracMatch",
    "fracGap",
    "score",
    "entroy",
    "pattern",
)
*/
pub fn read_trf_monomers(infile: impl AsRef<Path>) -> eyre::Result<MotifMonomers> {
    let reader = BufReader::new(File::open(infile)?);
    let mut motif_monomers: MotifMonomers = HashMap::new();
    for line in reader.lines().map_while(Result::ok) {
        let Some((
            motif,
            st,
            end,
            period,
            copy_num,
            _frac_match,
            _frac_gap,
            _score,
            _entropy,
            pattern,
        )) = line.split('\t').collect_tuple()
        else {
            continue;
        };

        let tname = motif.to_owned();
        let monomer = Monomer {
            srf_repeat: tname.clone(),
            trf_monomer: pattern.to_owned(),
            trf_period: period.parse()?,
            trf_copy_num: copy_num.parse()?,
        };
        let itv = Interval {
            start: st.parse()?,
            stop: end.parse()?,
            val: monomer,
        };
        motif_monomers
            .entry(tname.clone())
            .and_modify(|itree| itree.insert(itv.clone()))
            .or_insert_with(|| Lapper::new(vec![itv.clone()]));
    }
    Ok(motif_monomers)
}
