import gzip, sys
f_in, f_out, ncas, ncon = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
neff = round(4/(1/ncas+1/ncon)); kept=0
with gzip.open(f_in,'rt') as f, open(f_out,'w') as o:
    hdr=f.readline().rstrip('\n').split('\t'); i={c:k for k,c in enumerate(hdr)}
    o.write("SNP\tA1\tA2\tZ\tP\tN\n")
    for ln in f:
        p=ln.rstrip('\n').split('\t'); rs=p[i['rsids']]
        if not rs.startswith('rs'): continue
        try: b=float(p[i['beta']]); se=float(p[i['sebeta']]); pv=float(p[i['pval']])
        except: continue
        if se<=0 or not (0<pv<=1): continue
        o.write(f"{rs.split(',')[0]}\t{p[i['alt']].upper()}\t{p[i['ref']].upper()}\t{b/se:.5f}\t{pv:.6g}\t{neff}\n"); kept+=1
print(f"{f_out}: rows={kept:,} Neff={neff:,}")
