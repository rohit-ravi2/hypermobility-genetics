#!/usr/bin/env bash
# LD scores from the 1000G EUR MAGMA reference (MAF > 0.01, biallelic SNPs, 1 Mb window).
# Output goes to data/ref/ldsc_ld/, which run_ldsc.sh uses for both --ref-ld and --w-ld.
set -u
cd "$(dirname "$0")/.."
R=data/ref; O=$R/ldsc_ld; mkdir -p $O
for c in $(seq 1 22); do
  if [ ! -f $O/$c.l2.ldscore.gz ]; then
    "${PLINK2:-plink2}" --bfile $R/g1000_eur --chr $c --maf 0.01 --snps-only \
       --max-alleles 2 --make-bed --out $O/chr$c --silent 2>/dev/null
    python "${LDSC_DIR:-bin/ldsc}"/ldsc.py --bfile $O/chr$c --l2 --ld-wind-kb 1000 \
       --out $O/$c > $O/chr$c.ldsc.log 2>&1
    echo "$(date +%H:%M:%S) chr$c done: $(zcat $O/$c.l2.ldscore.gz 2>/dev/null | wc -l) snps"
    rm -f $O/chr$c.bed $O/chr$c.bim $O/chr$c.fam
  fi
done
echo LDSCORES_DONE
