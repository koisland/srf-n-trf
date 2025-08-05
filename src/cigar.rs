use eyre::{ContextCompat, bail};
use itertools::Itertools;
use paf::PafRecord;
use rust_lapper::Interval;

type IntervalPair = (Interval<u32, ()>, Interval<u32, ()>);

// \*|([0-9]+[MIDNSHP=X])+
#[derive(Debug, PartialEq, Eq)]
pub enum CigarToken {
    Number,
    Match,
    Mismatch,
    Insertion,
    Deletion,
    Softclip,
    Hardclip,
    Pad,
    Skip,
}

impl TryFrom<char> for CigarToken {
    type Error = eyre::Report;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        Ok(match value {
            '0'..='9' => CigarToken::Number,
            'M' => bail!("Ambiguous. Use extended cigar."),
            '=' => CigarToken::Match,
            'X' => CigarToken::Mismatch,
            'I' => CigarToken::Insertion,
            'D' => CigarToken::Deletion,
            'S' => CigarToken::Softclip,
            'H' => CigarToken::Hardclip,
            'P' => CigarToken::Pad,
            'N' => CigarToken::Skip,
            _ => bail!("Invalid token ({value})."),
        })
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),
    Mismatch(u32),
    Insertion(u32),
    Deletion(u32),
    Softclip(u32),
    Hardclip(u32),
    Pad(u32),
    Skip(u32),
}

pub fn parse_cigar(cg: &str) -> eyre::Result<Vec<CigarOp>> {
    let cg_iter = cg
        .trim_start_matches("cg:Z:")
        .chars()
        .chunk_by(|e| CigarToken::try_from(*e).unwrap());
    let mut cg_iter = cg_iter.into_iter().peekable();
    let mut cigar_ops = vec![];

    while let (Some((tk, mut elems)), Some((ntk, nelems))) = (cg_iter.next(), cg_iter.next()) {
        let num = elems.join("").parse::<u32>()?;
        let cg_op = match (&tk, &ntk) {
            (CigarToken::Number, CigarToken::Match) => CigarOp::Match(num),
            (CigarToken::Number, CigarToken::Mismatch) => CigarOp::Mismatch(num),
            (CigarToken::Number, CigarToken::Insertion) => CigarOp::Insertion(num),
            (CigarToken::Number, CigarToken::Deletion) => CigarOp::Deletion(num),
            (CigarToken::Number, CigarToken::Softclip) => CigarOp::Softclip(num),
            (CigarToken::Number, CigarToken::Hardclip) => CigarOp::Hardclip(num),
            (CigarToken::Number, CigarToken::Pad) => CigarOp::Pad(num),
            (CigarToken::Number, CigarToken::Skip) => CigarOp::Skip(num),
            _ => bail!(
                "Invalid cigar op ({tk:?}{:?}, {ntk:?}{:?})",
                elems.collect_vec(),
                nelems.collect_vec()
            ),
        };
        cigar_ops.push(cg_op);
    }
    Ok(cigar_ops)
}

/// Get intervals from query that align to target that meet some minimum length.
pub fn get_aligned_paired_itvs(
    rec: &PafRecord,
    min_length: u32,
) -> eyre::Result<Vec<IntervalPair>> {
    let mut pos: u32 = rec.target_start();
    let mut qpos: u32 = rec.query_start();
    let cg = rec.cg().context("Record has no cigar.")?;
    let cg_ops = parse_cigar(cg)?;

    let mut paired_itvs = vec![];
    for cg_op in cg_ops {
        let (q_adj, t_adj) = match cg_op {
            CigarOp::Match(l) | CigarOp::Mismatch(l) => (l, l),
            CigarOp::Insertion(l) | CigarOp::Softclip(l) => (l, 0),
            CigarOp::Deletion(l) => (0, l),
            CigarOp::Hardclip(_) => continue,
            CigarOp::Pad(l) => (l, 0),
            CigarOp::Skip(l) => (0, l),
        };
        if q_adj > min_length && t_adj > min_length {
            paired_itvs.push((
                Interval {
                    start: qpos,
                    stop: qpos + q_adj,
                    val: (),
                },
                Interval {
                    start: pos,
                    stop: pos + t_adj,
                    val: (),
                },
            ));
        }

        qpos += q_adj;
        pos += t_adj;
    }

    Ok(paired_itvs)
}
