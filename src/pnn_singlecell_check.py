import pandas as pd, numpy as np
a = pd.read_csv('data/expression/allen_m1_trimmed_means.csv', index_col=0)
print("Allen M1:", a.shape, "(genes x clusters)")
cols = list(a.columns)
pv   = [c for c in cols if 'PVALB' in c.upper()]
sst  = [c for c in cols if 'SST' in c.upper() and 'PVALB' not in c.upper()]
exc  = [c for c in cols if c.startswith('Exc')]
glia = [c for c in cols if any(k in c for k in ['Astro','Oligo','OPC','Micro','Endo','VLMC'])]
print(f"PVALB clusters={len(pv)} SST={len(sst)} Exc={len(exc)} glia/vascular={len(glia)}")

pnn = ['ACAN','BCAN','NCAN','VCAN','TNR','HAPLN1','HAPLN2','HAPLN3','HAPLN4','PTPRZ1']
bulk = pd.read_csv('results/pnn_partition.csv').set_index('symbol')

rows=[]
for g in pnn:
    if g not in a.index:
        rows.append((g,np.nan,np.nan,np.nan,np.nan,'NOT_IN_ALLEN',bulk.loc[g,'partition'] if g in bulk.index else '?')); continue
    r=a.loc[g]
    rows.append((g, r[pv].max(), r[pv].mean(), r[exc].mean(), r[glia].mean(),
                 'PV_ENRICHED' if r[pv].mean() > max(r[exc].mean(), r[glia].mean()) else 'not_PV',
                 bulk.loc[g,'partition'] if g in bulk.index else '?'))
df=pd.DataFrame(rows,columns=['gene','PV_max','PV_mean','Exc_mean','Glia_mean','sc_verdict','bulk_GTEx_call'])
pd.set_option('display.width',200)
print("\n=== PNN genes: Allen single-cell (trimmed means) vs bulk GTEx call ===")
print(df.to_string(index=False,float_format=lambda x:f"{x:,.2f}"))
df.to_csv('results/pnn_singlecell_vs_bulk.csv',index=False)

n_rescued = ((df['bulk_GTEx_call'].isin(['JOINT_ONLY','NEITHER'])) & (df['PV_max']>1)).sum()
print(f"\nGenes called non-brain by bulk GTEx but detectably expressed in PV interneurons: {n_rescued}")
