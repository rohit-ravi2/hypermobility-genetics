import numpy as np, pandas as pd

# gene locations (Entrez, chr, start, stop, strand, symbol)
gl = pd.read_csv('data/ref/NCBI38.gene.loc', sep='\t', header=None,
                 names=['entrez','chr','start','stop','strand','symbol'], dtype={'chr':str})
gl = gl[gl['chr'].isin([str(i) for i in range(1,23)])].copy()
gl['chr']=gl['chr'].astype(int)
gl['length']=gl['stop']-gl['start']

# SNP positions from the 1000G EUR reference MAGMA will actually use
bim = pd.read_csv('data/ref/g1000_eur.bim', sep='\t', header=None,
                  names=['chr','snp','cm','bp','a1','a2'], usecols=['chr','bp'],
                  dtype={'chr':str})
bim = bim[bim['chr'].isin([str(i) for i in range(1,23)])]
bim['chr']=bim['chr'].astype(int)
print(f"reference SNPs (autosomal): {len(bim):,}")

WIN=10000  # MAGMA-style window
snpcount=np.zeros(len(gl),dtype=int)
for c,sub in bim.groupby('chr'):
    pos=np.sort(sub['bp'].values)
    m=(gl['chr']==c).values
    idx=np.where(m)[0]
    lo=np.searchsorted(pos, gl['start'].values[m]-WIN,'left')
    hi=np.searchsorted(pos, gl['stop'].values[m]+WIN,'right')
    snpcount[idx]=hi-lo
gl['nsnp']=snpcount
gl['density']=gl['nsnp']/((gl['length']+2*WIN)/1000.0)   # SNPs per kb

def read_grp(n):
    return sorted({l.strip() for l in open(f'data/genesets/{n}.grp')
                   if l.strip() and not l.startswith('#') and l.strip()!=n})

sets={s:read_grp(s) for s in
      ['NABA_MATRISOME','NABA_CORE_MATRISOME','NABA_COLLAGENS','NABA_PROTEOGLYCANS',
       'NABA_ECM_GLYCOPROTEINS','GOCC_PERINEURONAL_NET']}

part=pd.read_csv('results/core_matrisome_partition.csv')
for p in ['SHARED','JOINT_ONLY','BRAIN_ONLY']:
    sets[f'CORE_{p}']=sorted(part.loc[part['partition']==p,'input'])

bg_len=gl['length'].median(); bg_n=gl['nsnp'].median(); bg_d=gl['density'].median()
print(f"\nGENOME BACKGROUND (n={len(gl):,} genes): median length={bg_len:,.0f} bp  "
      f"median nSNP={bg_n:.0f}  median density={bg_d:.2f}/kb\n")

rows=[]
for name,genes in sets.items():
    sub=gl[gl['symbol'].isin(genes)]
    if len(sub)<3: continue
    rows.append(dict(set=name, n_matched=len(sub), n_input=len(genes),
        med_len=sub['length'].median(), len_ratio=sub['length'].median()/bg_len,
        med_nsnp=sub['nsnp'].median(), nsnp_ratio=sub['nsnp'].median()/bg_n,
        med_density=sub['density'].median(), dens_ratio=sub['density'].median()/bg_d))
r=pd.DataFrame(rows).sort_values('len_ratio',ascending=False)
pd.set_option('display.width',220)
print(r.to_string(index=False,float_format=lambda x:f"{x:,.2f}"))
r.to_csv('results/gene_length_density_diagnostic.csv',index=False)

from scipy import stats
print("\n-- Mann-Whitney U: set gene length vs all other genes --")
for name,genes in sets.items():
    sub=gl[gl['symbol'].isin(genes)]
    oth=gl[~gl['symbol'].isin(genes)]
    if len(sub)<3: continue
    u,p=stats.mannwhitneyu(sub['length'],oth['length'],alternative='greater')
    print(f"{name:28s} n={len(sub):4d}  median={sub['length'].median():>10,.0f}  p(longer)={p:.3g}")
