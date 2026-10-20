import pandas as pd, numpy as np
gl=pd.read_csv('data/ref/NCBI38.gene.loc',sep='\t',header=None,
               names=['entrez','chr','start','stop','strand','symbol'],dtype=str)
sym2ez=dict(zip(gl['symbol'],gl['entrez']))

def grp(n):
    return sorted({l.strip() for l in open(f'data/genesets/{n}.grp')
                   if l.strip() and not l.startswith('#') and l.strip()!=n})

part=pd.read_csv('results/core_matrisome_partition.csv')
sets={}
for p in ['SHARED','JOINT_ONLY','BRAIN_ONLY']:
    sets[f'CORE_{p}']=sorted(part.loc[part['partition']==p,'input'])
sets['PV_MATRISOME']=grp('PV_MATRISOME')
sets['PNN_GO']=grp('GOCC_PERINEURONAL_NET')
for s in ['NABA_CORE_MATRISOME','NABA_MATRISOME','NABA_COLLAGENS','NABA_PROTEOGLYCANS']:
    sets[s]=grp(s)

# random length-matched null sets drawn from all genes
gl2=gl.copy(); gl2['len']=gl2['stop'].astype(int)-gl2['start'].astype(int)
rng=np.random.default_rng(20260814)
tgt=gl2[gl2['symbol'].isin(sets['PV_MATRISOME'])]['len'].median()
pool=gl2[(gl2['len']>tgt*0.5)&(gl2['len']<tgt*2.0)]
for i in range(5):
    sets[f'RANDOM_LENMATCH_{i+1}']=list(pool.sample(len(sets['PV_MATRISOME']),random_state=100+i)['symbol'])

with open('genesets/genesets.setannot','w') as f:
    for name,genes in sets.items():
        ez=[sym2ez[g] for g in genes if g in sym2ez]
        print(f"{name:26s} input={len(genes):5d} mapped={len(ez):5d}")
        f.write(name+' '+' '.join(ez)+'\n')
