import gzip, re, sys
import pandas as pd, numpy as np

GS='data/genesets'; EX='data/expression'

def read_grp(name):
    genes=[]
    for ln in open(f'{GS}/{name}.grp'):
        ln=ln.strip()
        if not ln or ln.startswith('#') or ln==name: continue
        genes.append(ln)
    return sorted(set(genes))

pnn_core = read_grp('GOCC_PERINEURONAL_NET')
core     = read_grp('NABA_CORE_MATRISOME')
matri    = read_grp('NABA_MATRISOME')
print(f"PNN_CORE({len(pnn_core)}): {' '.join(pnn_core)}")
print(f"NABA_CORE_MATRISOME: {len(core)}   NABA_MATRISOME: {len(matri)}")

# ---- HGNC symbol -> entrez, plus alias resolution ----
h = pd.read_csv(f'{GS}/hgnc_complete_set.txt', sep='\t', dtype=str, low_memory=False)
h = h[h['status']=='Approved']
sym2ez = dict(zip(h['symbol'], h['entrez_id']))
alias2sym={}
for col in ['alias_symbol','prev_symbol']:
    for s,a in zip(h['symbol'], h[col].fillna('')):
        for x in a.split('|'):
            if x: alias2sym.setdefault(x, s)
def canon(g):
    if g in sym2ez: return g
    return alias2sym.get(g, g)

# ---- GTEx v11 ----
gt = pd.read_csv(f'{EX}/gtex_v11_median_tpm.gct.gz', sep='\t', skiprows=2, low_memory=False)
gt['Description']=gt['Description'].astype(str)
brain = [c for c in gt.columns if c.startswith('Brain_')]
conn  = ['Cells_Cultured_fibroblasts','Skin_Not_Sun_Exposed_Suprapubic','Skin_Sun_Exposed_Lower_leg',
         'Artery_Aorta','Artery_Coronary','Artery_Tibial','Adipose_Subcutaneous','Adipose_Visceral_Omentum']
conn  = [c for c in conn if c in gt.columns]
print(f"\nGTEx: {len(brain)} brain cols, {len(conn)} connective cols")

gt['brain_max']=gt[brain].max(axis=1); gt['conn_max']=gt[conn].max(axis=1)
expr = gt.groupby('Description')[['brain_max','conn_max']].max()

TH=1.0
def partition(genes, label):
    rows=[]
    for g in genes:
        cg=canon(g)
        key = cg if cg in expr.index else (g if g in expr.index else None)
        if key is None:
            rows.append((g,cg,np.nan,np.nan,'NO_EXPR_DATA')); continue
        b,c = expr.loc[key,'brain_max'], expr.loc[key,'conn_max']
        if b>=TH and c>=TH: p='SHARED'
        elif b>=TH:         p='BRAIN_ONLY'
        elif c>=TH:         p='JOINT_ONLY'
        else:               p='NEITHER'
        rows.append((g,cg,b,c,p))
    df=pd.DataFrame(rows,columns=['input','symbol','brain_max','conn_max','partition'])
    print(f"\n=== {label} (n={len(df)}) ===")
    print(df['partition'].value_counts().to_string())
    return df

pnn_df  = partition(pnn_core,'PNN_CORE')
print(pnn_df[['symbol','brain_max','conn_max','partition']].to_string(index=False))
core_df = partition(core,'NABA_CORE_MATRISOME')
core_df.to_csv('results/core_matrisome_partition.csv',index=False)
pnn_df.to_csv('results/pnn_partition.csv',index=False)

# EDS monogenic genes - the GAG vs collagen split
eds = {'COL5A1':'cEDS','COL5A2':'cEDS','COL1A1':'aEDS','COL1A2':'cvEDS','COL3A1':'vEDS',
       'TNXB':'clEDS','AEBP1':'clEDS2','ADAMTS2':'dEDS','PLOD1':'kEDS','FKBP14':'kEDS',
       'ZNF469':'BCS','PRDM5':'BCS','B4GALT7':'spEDS','B3GALT6':'spEDS','SLC39A13':'spEDS',
       'CHST14':'mcEDS','DSE':'mcEDS','COL12A1':'mEDS','C1R':'pEDS','C1S':'pEDS',
       'FBN1':'Marfan','FBN2':'CCA'}
rows=[]
for g,sub in eds.items():
    cg=canon(g)
    if cg in expr.index:
        rows.append((g,sub,expr.loc[cg,'brain_max'],expr.loc[cg,'conn_max']))
e=pd.DataFrame(rows,columns=['gene','subtype','brain_max','conn_max']).sort_values('brain_max',ascending=False)
e['brain_conn_ratio']=(e['brain_max']/e['conn_max'].replace(0,np.nan)).round(3)
print("\n=== EDS/monogenic connective-tissue genes: GTEx v11 brain vs connective (max TPM) ===")
print(e.to_string(index=False))
e.to_csv('results/eds_gene_expression.csv',index=False)
