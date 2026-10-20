#!/usr/bin/env bash
# Per-trait MAGMA gene analysis (10 kb windows, 1000G EUR). Expects work/pval/<TRAIT>.tsv
# (SNP, P, N) and work/annot_10kb.genes.annot from `magma --annotate window=10,10`.
cd "$(dirname "$0")/.."
for T in ADHD ASD ANX MDD; do
 ( "${MAGMA:-magma}" --bfile data/ref/g1000_eur \
     --pval work/pval/$T.tsv use=SNP,P ncol=N \
     --gene-annot work/annot_10kb.genes.annot \
     --out results/magma_$T > work/magma_$T.log 2>&1 ; echo "$T exit=$?" ) &
done
wait
echo ALL_GENE_ANALYSES_DONE

# Competitive gene-set analysis on the sets in genesets/genesets.setannot
# (written by src/make_setfile.py).
for T in ADHD ASD ANX MDD; do
  "${MAGMA:-magma}" --gene-results results/magma_$T.genes.raw \
    --set-annot genesets/genesets.setannot --out results/gsa_$T
done
