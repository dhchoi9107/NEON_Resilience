# Recompute rarefied LCBD for ALL 26 sites at RARE_N=10, N_DRAWS=999 (new standard).
# Validated against 30-draw pipeline. Updates LCBD_*_rare cols in FINAL_v2_pooled_26 & lidar_pooled_predictors_26.
import os, shutil, numpy as np, pandas as pd
VS="E:/neon_lidar/vegetation_structure"; D="NEON_v2/data"
RARE_N=10; NDRAWS=999; SEED=11
mt=pd.read_csv(f"{VS}/vst_mappingandtagging.csv",usecols=['individualID','plotID','taxonID','siteID'],low_memory=False)
ai=pd.read_csv(f"{VS}/vst_apparentindividual.csv",usecols=['individualID','plantStatus','stemDiameter'],low_memory=False)
df=ai.merge(mt,on='individualID',how='left')
df=df[df.plantStatus.astype(str).str.contains('Live',na=False)]
df=df[df.stemDiameter>=10]
df=df[df.taxonID.notna() & ~df.taxonID.astype(str).str.contains('2PLANT|UNK',na=False)].drop_duplicates('individualID')
def hill_cov(c):
    c=np.asarray([x for x in c if x>0]); N=c.sum()
    if N==0: return 0,np.nan
    f1=int((c==1).sum()); f2=int((c==2).sum())
    return int(N),(1-(f1/N)*(((N-1)*f1)/((N-1)*f1+2*max(f2,1))) if N>1 else np.nan)
cr=[]
for pid,g in df.groupby('plotID'):
    N,cov=hill_cov(g.groupby('taxonID').size().values); cr.append((pid,g.siteID.iloc[0],N,cov))
cv=pd.DataFrame(cr,columns=['plotID','siteID','N','cvg'])
def baselga_vec(pa):
    pa=pa.astype(float); a=pa@pa.T; notpa=1-pa; b=pa@notpa.T; c=notpa@pa.T
    mbc=np.minimum(b,c); den=a+mbc
    sim=np.where(den>0,mbc/den,0.0); sor=np.where((2*a+b+c)>0,(b+c)/(2*a+b+c),0.0)
    return sim,sor-sim
good=set(cv[(cv['cvg']>=0.9)&(cv['N']>=RARE_N)]['plotID'])
rng=np.random.RandomState(SEED); rows=[]
for site,g in df[df.plotID.isin(good)].groupby('siteID'):
    piv=g.groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    if piv.shape[0]<2: continue
    plots=piv.index.tolist(); base=piv.values.astype(int); n=len(plots); Sn=base.shape[1]
    pools=[np.repeat(np.arange(Sn),base[i]) for i in range(n)]; tu=np.zeros(n); ne=np.zeros(n)
    for _ in range(NDRAWS):
        rar=np.zeros((n,Sn))
        for i in range(n): rar[i]=np.bincount(rng.choice(pools[i],RARE_N,replace=False),minlength=Sn)
        sim,sne=baselga_vec(rar>0); off=~np.eye(n,dtype=bool)
        tu+=(sim*off).sum(1)/(n-1); ne+=(sne*off).sum(1)/(n-1)
    for i,p in enumerate(plots): rows.append((p,tu[i]/NDRAWS,ne[i]/NDRAWS))
lc=pd.DataFrame(rows,columns=['plotID','LCBD_turnover_rare','LCBD_nestedness_rare'])
print(f"999-draw rare LCBD computed: {len(lc)} plots")
# update both files
for fn in ["FINAL_v2_pooled_26.csv","lidar_pooled_predictors_26.csv"]:
    p=os.path.join(D,fn); d=pd.read_csv(p)
    if 'LCBD_turnover_rare' not in d.columns: print(f"  {fn}: no rare cols, skip"); continue
    shutil.copy2(p,p+".bak3")
    old=d[['plotID','LCBD_turnover_rare','LCBD_nestedness_rare']].copy()
    d=d.drop(columns=['LCBD_turnover_rare','LCBD_nestedness_rare']).merge(lc,on='plotID',how='left')
    # keep original column order
    d.to_csv(p,index=False)
    m=old.merge(lc,on='plotID',suffixes=('_old','_new')).dropna()
    print(f"  {fn}: updated | turnover Δmedian {(m.LCBD_turnover_rare_new-m.LCBD_turnover_rare_old).abs().median():.4f}, nestedness Δmedian {(m.LCBD_nestedness_rare_new-m.LCBD_nestedness_rare_old).abs().median():.4f}")
print("DONE")
