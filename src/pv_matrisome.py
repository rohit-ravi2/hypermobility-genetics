import pandas as pd, numpy as np
a = pd.read_csv('data/expression/allen_m1_trimmed_means.csv', index_col=0)
cols=list(a.columns)
pv   = [c for c in cols if 'PVALB' in c.upper()]
other= [c for c in cols if c not in pv]

def read_grp(n):
    return sorted({l.strip() for l in open(f'data/genesets/{n}.grp')
                   if l.strip() and not l.startswith('#') and l.strip()!=n})
matri=read_grp('NABA_MATRISOME'); core=read_grp('NABA_CORE_MATRISOME')

present=[g for g in matri if g in a.index]
print(f"matrisome genes in Allen M1: {len(present)}/{len(matri)}")
sub=a.loc[present]
pv_mean=sub[pv].mean(axis=1); oth_max=sub[other].max(axis=1); oth_mean=sub[other].mean(axis=1)
res=pd.DataFrame({'PV_mean':pv_mean,'PV_max':sub[pv].max(axis=1),
                  'other_mean':oth_mean,'other_max':oth_max})
res['spec']=res['PV_mean']/(res['other_mean']+0.05)
res['is_core']=[g in core for g in res.index]

# expressed in PV at all, and enriched vs the rest of cortex
sel=res[(res['PV_mean']>=0.25)&(res['spec']>=1.5)].sort_values('spec',ascending=False)
print(f"\n=== PV-enriched matrisome genes (PV_mean>=0.25, specificity>=1.5): n={len(sel)} ===")
pd.set_option('display.width',200)
print(sel.head(30).to_string(float_format=lambda x:f"{x:,.2f}"))
sel.to_csv('results/pv_matrisome_geneset.csv')
open('data/genesets/PV_MATRISOME.grp','w').write('\n'.join(sel.index)+'\n')
print(f"\ncore-matrisome members among them: {sel['is_core'].sum()}")
print("PNN canon present:", [g for g in ['ACAN','HAPLN1','HAPLN4','BCAN','TNR','NCAN','VCAN','PTPRZ1'] if g in sel.index])
