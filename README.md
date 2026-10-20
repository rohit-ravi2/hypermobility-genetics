# hypermobility-genetics

Common-variant genetics of joint hypermobility from public GWAS. Using FinnGen R13 (ICD-10 M35.7,
2,199 cases / 483,310 controls) I estimate the SNP heritability of clinically diagnosed
hypermobility at h2 = 0.088 (SE 0.018, liability scale) and find positive genetic correlations
with hernia (rg = 0.24), varicose veins (0.28) and pelvic organ prolapse (0.26). A PIEZO1 variant
previously reported for Beighton score replicates in FinnGen (p = 0.0026). The other half of the
project is a negative result: extracellular-matrix and perineuronal-net gene sets carry no
excess common-variant signal for ADHD, autism, anxiety or depression.

Write-up: https://www.rohitravi.com/projects/hypermobility-brain

## Results

### Heritability (LDSC)

`results/ldsc/h2_hyper_liab.log`, `results/ldsc/h2_hypermobility.log`

| | estimate |
|---|---|
| h2, liability scale (population prevalence 0.45%) | 0.088 (0.0176) |
| h2, observed scale | 0.19 (0.038) |
| intercept | 1.025 (0.004) |
| ratio | 0.44 (0.07) |
| mean chi2 | 1.058 |

LD scores were computed from the 1000G EUR panel (`scripts/build_ldscores.sh`). As a check on
them, common executive function (Hatoum et al. 2023) gives h2 = 0.0947 (SE 0.0063) through the
same pipeline (`results/ldsc/h2_execfunction.log`).

### Genetic correlations

`results/ldsc/rg_hyper_panel.log`, HapMap3 SNPs, FinnGen R13 endpoints:

| hypermobility vs | rg | SE | p |
|---|---|---|---|
| varicose veins (I9_VARICVE) | 0.278 | 0.087 | 0.0015 |
| female genital prolapse (N14_FEMGENPROL) | 0.259 | 0.092 | 0.0048 |
| hernia (K11_HERNIA) | 0.240 | 0.094 | 0.011 |
| common executive function (Hatoum 2023) | -0.124 | 0.105 | 0.24 |

Among the three surgical traits (`results/ldsc/rg_lax_*.log`): hernia-prolapse rg = 0.437
(0.043), hernia-varicose 0.148 (0.041), prolapse-varicose 0.114 (0.037). Varicose veins share
little with the other two, but hypermobility correlates with all three at a similar level.

Cross-trait intercepts for the hypermobility pairs range from -0.013 to 0.033.

### PIEZO1 rs78579285

Reported for Beighton score at p = 6e-12 in a single 23andMe scan with no public summary
statistics. In FinnGen R13 hypermobility it gives beta = -0.120 (SE 0.040), p = 0.0026, a
nominal replication of one pre-specified variant. The same variant in FinnGen varicose veins:
beta = -0.125, p = 1.7e-33. `results/piezo1_rs78579285_finngen_r13.tsv`.

### Matrisome and perineuronal-net gene sets (MAGMA)

The hypothesis was that connective-tissue genes and perineuronal-net genes share a cause, so
matrisome/PNN sets should carry excess common-variant signal for neurodevelopmental and
psychiatric traits. MAGMA 1.10 competitive gene-set analysis, conditioned on gene size and
density, on ADHD (Demontis 2023), autism (Grove 2019), anxiety (PGC 2026) and depression
(PGC 2025). Sets: NABA matrisome divisions, GO perineuronal net, the core matrisome split by
GTEx expression into brain-only / joint-only / shared, a PV-interneuron matrisome set from Allen
M1 single-cell data, and five random sets matched on gene length.

- Nothing survives correction (36 primary tests). The smallest p is 0.034 (PNN_GO, depression,
  10 genes).
- A random length-matched set beat the best candidate on the same trait (autism: random
  p = 0.047, PV_MATRISOME p = 0.050).
- NABA_MATRISOME (970-1,003 genes) is null on all four traits, p = 0.77-0.97, with negative betas.
  Earlier claims of hundreds of matrisome genes associated with autism came from catalogue
  counting without a gene-size correction.
- Collagen genes have a median length 5.5x the genome median (`results/gene_length_density_diagnostic.csv`),
  so uncorrected gene-based tests on these sets inflate.

Results: `results/magma/gsa_<TRAIT>.gsa.out`.

Supporting tables:
- `results/eds_gene_expression.csv`: GTEx v11 brain vs connective-tissue expression of the
  monogenic EDS genes. The spondylodysplastic EDS genes (SLC39A13, B4GALT7, B3GALT6) have
  brain/connective expression ratios of 0.50-0.75; the collagen genes are at 0.001-0.04.
- `results/pnn_singlecell_vs_bulk.csv`: the GO perineuronal-net set checked against PV
  interneurons in Allen M1 single-cell data. Only 3 of its 10 genes (ACAN, HAPLN1, HAPLN4) are
  PV-enriched.

## Repo layout

```
src/
  build_genesets.py             partition core matrisome by GTEx brain vs connective expression
  pv_matrisome.py               PV-interneuron-enriched matrisome set (Allen M1)
  pnn_singlecell_check.py       is the GO PNN set actually PV-specific
  length_density_diagnostic.py  gene length / SNP density of each set vs genome
  make_setfile.py               writes genesets/genesets.setannot for MAGMA
  prep_finngen.py               FinnGen sumstats -> LDSC input with effective N
scripts/
  fetch_sumstats.sh             downloads the GWAS
  build_ldscores.sh             1000G EUR LD scores
  run_ldsc.sh                   munge, h2 and rg calls
  run_magma_gene.sh             MAGMA gene and gene-set analysis
genesets/                       MAGMA set file and the derived PV_MATRISOME list
results/                        tables; ldsc/ and magma/ hold the run logs and outputs
```

## How to run

Needs LDSC (v3.0.1 was used), MAGMA 1.10, PLINK 2 and Python with pandas,
numpy and scipy. Run from the repo root.

```
bash scripts/fetch_sumstats.sh
bash scripts/build_ldscores.sh
bash scripts/run_ldsc.sh
python src/build_genesets.py && python src/pv_matrisome.py && python src/make_setfile.py
bash scripts/run_magma_gene.sh
```

`run_magma_gene.sh` expects per-trait SNP/P/N files in `work/pval/` and the 10 kb annotation in
`work/annot_10kb.genes.annot`; the column layout differs for each PGC file, so I did that step by
hand.

## Data

Not included. All public:

- FinnGen R13 summary statistics (M13_HYPERMOBILITY, K11_HERNIA, I9_VARICVE, N14_FEMGENPROL):
  https://www.finngen.fi/en/access_results
- PGC: ADHD (Demontis 2023), autism (Grove 2019), anxiety (2026), depression (2025, no 23andMe),
  https://pgc.unc.edu/for-researchers/download-results/
- Common executive function, Hatoum et al. 2023, GWAS Catalog GCST90162547.
- 1000G EUR reference and NCBI gene locations from the MAGMA site; HapMap3 SNP list from the
  LDSC resources (saved as `data/ref/merge_alleles.txt`).
- Gene sets: MSigDB (NABA_*, GOCC_PERINEURONAL_NET), HGNC complete set. Expression: GTEx v11
  median TPM, Allen Brain Map human M1 trimmed means.

`scripts/fetch_sumstats.sh` has the exact URLs for the GWAS files.

## Notes

- The MAGMA SNP annotation used SNP positions from the 1000G EUR reference (GRCh37) with
  `NCBI38.gene.loc` (GRCh38) gene boundaries, so some SNP-to-gene assignments are off. The
  gene-set results above should be re-run with `NCBI37.3.gene.loc` before being relied on. The
  LDSC results match on rsID and are not affected.
- Hypermobility here is an ICD-10 registry diagnosis, not a Beighton score or hEDS by the 2017
  criteria, and the exposure GWAS is small (N_eff 8,756). Its h2 z-score is about 5.

## License

MIT
