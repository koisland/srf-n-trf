# `srf-n-trf`
Script to take `srf` and `trf` output and produce a BED file with only regions corresponding to monomers of a given periodicity.

## Getting Started
1. Run `srf` and `trf` on a T2T genome assembly using https://github.com/logsdon-lab/Snakemake-srf.git. Use these outputs.
```
results/{sample}/{contig}
├── count.txt
├── monomers.tsv (*)
├── srf.bed
├── srf.fa
└── srf.paf (*)
```

2. Clone repo. Requires `rust`.
```bash
git clone https://github.com/koisland/srf-n-trf.git
```

3. Compile.
```bash
cargo build --release
```

## Usage
```bash
srf-n-trf \
-p results/{sample}/{contig}/srf.paf \
-m results/{sample}/{contig}/monomers.tsv \
-s 170 340 42 \
-d 0.02
```
This will:
* Search for `trf` monomers in PAF cigar with a periodicity of `170`, `340`, and `42` with a `2%` length difference.
    * These correspond to α-satellite and HSAT-1A repeats.
* Generate a BED9 file in target coordinate space and the overlapping monomers delimited by commas in the `name` column.

```
chr3_mat_hsa4	76301546	76301589	TATGAAAAGAAAGGTTAAACTCTGTGAGTTGAACGCACACATCACAAAGTAGTTTCTGAGAATGATTCTCTCTAGTTTTTATACGAAGATATTTCCTTTTCTACCATTGGCCTCAAAGCACTTGAAATCTCCACCTGCAAATTCCACAAAAAGAGTGTTTCAAATCTGCTCTGTCTAAAGGAAGCTTCAACTCTGTGAGTTGAATACACACAACACAAAGAAGTTACTGAGAATTCTTCTGTCTAGCATTATATGAAGAAATCCCGTTTCCAACGAAGGCCTCAAAGAGGTCCAAATATCCACTTGCAGACTTAACAAACAGAGTGTTTCCAAACTGCTC,AAAAGAAAGGTTAAACTCTGTGAGTTGAACACACACAACACAAAGAAGTTACTGAGAATGATTCTGTCTAGCATTATACGAAGAAATCCCGTTTCCAACGAAGGCCTCAAAGAGGTCCAAATATCCACTTGCAACTTAACAAACAGAGTGTTTCCAAACTGCTCTGTC,AAAGGAAGGTTCAACTCTGTGAGTTGAACACACACATCACAAAGAAGTTACTGAGAATGATTCTCTCTAGTTTTATACGAAGATATTTCCTTTTCAAAAATGGCCTCAAAGCGCTTCAAATCTCCACTTGCAAATTCCACAAAAAGAGTGTTTCAAATCTGCTCTGTCT	0	-	76301546	76301589	0,0,0
```
> The example in `test/chr3` is mGorGor1 chr3_mat_hsa4 from the [T2T Primates project](https://github.com/marbl/Primates?tab=readme-ov-file).

## Examples
```bash
cargo run -- -p <(zcat test/chr3/srf.paf.gz) -m <(zcat test/chr3/monomers.tsv) # mGorGor1
cargo run -- -p <(zcat test/chrX/srf.paf.gz) -m <(zcat test/chrX/monomers.tsv) # mPonAbe1
```

## TODO
* [ ] Unit and integration tests.
* [ ] Support compressed output.
* [ ] Better documentation.
* [ ] CI and release workflow.
