#!/usr/bin/env bash
# Munging, h2 and rg calls behind the logs in results/ldsc/ and results/munge_logs/.
# Needs LDSC (set LDSC_DIR), data/ref/ldsc_ld/ from build_ldscores.sh and
# data/ref/merge_alleles.txt (HapMap3 w_hm3.snplist).
set -eu
cd "$(dirname "$0")/.."
L=${LDSC_DIR:-bin/ldsc}
mkdir -p work results/ldsc

# FinnGen -> LDSC input with effective N = 4 / (1/cases + 1/controls).
# Hypermobility: 2,199 cases / 483,310 controls -> Neff 8,756. Neff used for the
# other traits (from the FinnGen R13 case/control counts): hernia 260,115,
# varicose veins 143,026, female genital prolapse 80,454.
python src/prep_finngen.py data/sumstats/finngen_R13_hypermobility.gz work/hypermobility.forldsc.tsv 2199 483310
# python src/prep_finngen.py data/sumstats/finngen_R13_<trait>.gz work/<trait>.forldsc.tsv <cases> <controls>
# work/execfunction.forldsc.tsv is the Hatoum 2023 file reduced to SNP A1 A2 Z P N (Z = beta/se).

for T in hypermobility hernia varicose prolapse execfunction; do
  python $L/munge_sumstats.py --sumstats work/$T.forldsc.tsv --out work/${T}_hm3 \
    --merge-alleles data/ref/merge_alleles.txt --chunksize 500000 --signed-sumstats Z,0
done
for T in hypermobility execfunction; do
  python $L/munge_sumstats.py --sumstats work/$T.forldsc.tsv --out work/$T \
    --chunksize 500000 --signed-sumstats Z,0
done

LD="--ref-ld-chr data/ref/ldsc_ld/ --w-ld-chr data/ref/ldsc_ld/"
python $L/ldsc.py --out results/ldsc/h2_hypermobility --h2 work/hypermobility.sumstats.gz $LD
python $L/ldsc.py --out results/ldsc/h2_hyper_liab --h2 work/hypermobility.sumstats.gz $LD \
  --samp-prev 0.5 --pop-prev 0.0045
python $L/ldsc.py --out results/ldsc/h2_execfunction --h2 work/execfunction.sumstats.gz $LD

python $L/ldsc.py --out results/ldsc/rg_hyper_panel $LD --rg \
  work/hypermobility_hm3.sumstats.gz,work/hernia_hm3.sumstats.gz,work/varicose_hm3.sumstats.gz,work/prolapse_hm3.sumstats.gz,work/execfunction_hm3.sumstats.gz
for T in hernia varicose prolapse; do
  python $L/ldsc.py --out results/ldsc/rg_lax_$T $LD --rg \
    work/${T}_hm3.sumstats.gz,work/hernia_hm3.sumstats.gz,work/varicose_hm3.sumstats.gz,work/prolapse_hm3.sumstats.gz
done
