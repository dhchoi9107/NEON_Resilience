import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
VS="E:/neon_lidar/vegetation_structure"; D="NEON_v2/data"
NDRAWS=999; SEED=11
# --- pooled SAD, all 26 sites (same filters as build_new_pooled) ---
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
    cov=1-(f1/N)*(((N-1)*f1)/((N-1)*f1+2*max(f2,1))) if N>1 else np.nan
    return int(N),cov
# per-plot N & coverage
covrows=[]
for pid,g in df.groupby('plotID'):
    N,cov=hill_cov(g.groupby('taxonID').size().values)
    covrows.append((pid,pid[:4] if False else g.siteID.iloc[0],N,cov))
cv=pd.DataFrame(covrows,columns=['plotID','siteID','N','cov'])

def baselga_vec(pa):  # pa: n x S presence-absence -> sim(turnover), sne(nestedness) matrices
    pa=pa.astype(float); a=pa@pa.T; notpa=1-pa
    b=pa@notpa.T; c=notpa@pa.T
    mbc=np.minimum(b,c); denom=a+mbc
    sim=np.where(denom>0,mbc/denom,0.0)
    sor=np.where((2*a+b+c)>0,(b+c)/(2*a+b+c),0.0)
    sne=sor-sim
    return sim,sne

def rare_lcbd(RARE_N):
    good=set(cv[(cv['cov']>=0.9)&(cv['N']>=RARE_N)]['plotID'])
    rng=np.random.RandomState(SEED); rows=[]
    for site,g in df[df.plotID.isin(good)].groupby('siteID'):
        piv=g.groupby(['plotID','taxonID']).size().unstack(fill_value=0)
        if piv.shape[0]<2: continue
        plots=piv.index.tolist(); base=piv.values.astype(int); n=len(plots); S=base.shape[1]
        pools=[np.repeat(np.arange(S),base[i]) for i in range(n)]
        tu=np.zeros(n); ne=np.zeros(n)
        for _ in range(NDRAWS):
            rar=np.zeros((n,S))
            for i in range(n):
                pick=rng.choice(pools[i],RARE_N,replace=False)
                rar[i]=np.bincount(pick,minlength=S)
            sim,sne=baselga_vec(rar>0)
            off=~np.eye(n,dtype=bool)
            tu+=(sim*off).sum(1)/(n-1); ne+=(sne*off).sum(1)/(n-1)
        for i,p in enumerate(plots):
            rows.append((p,tu[i]/NDRAWS,ne[i]/NDRAWS))
    return pd.DataFrame(rows,columns=['plotID','LCBD_turnover_rare','LCBD_nestedness_rare'])

# --- 157 model for beta responses ---
base=pd.read_csv(f"{D}/FINAL_v2_pooled_26.csv")
for f in ["plot_modis_gpp_26.csv","plot_pml_gpp_26.csv"]: base=base.merge(pd.read_csv(f"{D}/{f}"),on="plotID",how="left")
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];G=["modis_gpp","pml_gpp"]
keep=S+P+Dn+G
def run157(base2,resp):
    d=base2[[resp,"domain","siteID"]+keep].dropna().copy()
    for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
    y=d[resp].astype(float); d["y"]=(y-y.mean())/y.std()
    f="y ~ C(domain)"+"".join(" + "+c for c in keep)
    m=smf.ols(f,d).fit(cov_type="cluster",cov_kwds={"groups":d["siteID"]})
    wd=m.wald_test(",".join(f"{t}=0" for t in Dn),use_f=True,scalar=True)
    wf=m.wald_test(",".join(f"{t}=0" for t in G),use_f=True,scalar=True)
    o=lambda pr: smf.ols("y ~ C(domain)"+"".join(" + "+c for c in pr),d).fit().rsquared
    return dict(n=len(d),dyn_p=float(wd.pvalue),dyn_dR2=o(S+P+Dn)-o(S+P),
                fun_p=float(wf.pvalue),fun_dR2=o(S+P+Dn+G)-o(S+P+Dn))

print(f"{'RARE_N':>7} {'resp':>11} {'n':>4} {'dyn_dR2':>8} {'dyn_p':>8} {'fun_dR2':>8} {'fun_p':>8}")
for RN in [10,15,20]:
    lc=rare_lcbd(RN)
    b2=base.drop(columns=['LCBD_turnover_rare','LCBD_nestedness_rare']).merge(lc,on='plotID',how='left')
    b2=b2[b2.sample_coverage>=0.9]
    for resp in ['LCBD_turnover_rare','LCBD_nestedness_rare']:
        r=run157(b2,resp)
        print(f"{RN:7d} {resp.replace('LCBD_','').replace('_rare',''):>11} {r['n']:4d} {r['dyn_dR2']:8.4f} {r['dyn_p']:8.4f} {r['fun_dR2']:8.4f} {r['fun_p']:8.4f}")
