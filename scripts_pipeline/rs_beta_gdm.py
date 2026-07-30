# "Saturation model" for compositional dissimilarity (Bray-Curtis in [0,1]) ~ structural + spectral distance.
# Beta regression (logit link, handles the [0,1] ceiling) with within-site restricted-permutation p (999).
import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.spatial.distance import pdist, squareform
import statsmodels.api as sm
VS="E:/neon_lidar/vegetation_structure"; D="NEON_v2/data/"
mt=pd.read_csv(f"{VS}/vst_mappingandtagging.csv",usecols=['individualID','plotID','taxonID','siteID'],low_memory=False)
ai=pd.read_csv(f"{VS}/vst_apparentindividual.csv",usecols=['individualID','plantStatus','stemDiameter'],low_memory=False)
df=ai.merge(mt,on='individualID',how='left')
df=df[df.plantStatus.astype(str).str.contains('Live',na=False)]; df=df[df.stemDiameter>=10]
df=df[df.taxonID.notna() & ~df.taxonID.astype(str).str.contains('2PLANT|UNK',na=False)].drop_duplicates('individualID')
f=pd.read_csv(D+"FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()
Scols=[c for c in ['VCI_mean','LAI_mean','Rugosity_mean','Vert_CV_mean','Canopy_Ht_mean','FHD_mean','Deep_Gap_mean'] if c in f.columns]
Pcols=[c for c in ['NDVI_mean','EVI_mean','ARVI_mean','SAVI_mean'] if c in f.columns]
for c in Scols+Pcols: f['z_'+c]=(f[c]-f[c].mean())/f[c].std()
sites=[]
for s,g in f.groupby('siteID'):
    plots=g.plotID.tolist(); sad=df[(df.siteID==s)&(df.plotID.isin(plots))].groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    common=[p for p in plots if p in sad.index]; gg=g.set_index('plotID').loc[common]
    if len(common)<6 or gg[['z_'+c for c in Scols]].isna().any().any() or gg[['z_'+c for c in Pcols]].isna().any().any(): continue
    sites.append((squareform(pdist(sad.loc[common].values,metric='braycurtis')),squareform(pdist(gg[['z_'+c for c in Scols]].values)),squareform(pdist(gg[['z_'+c for c in Pcols]].values)),len(common)))
iu=lambda M: M[np.triu_indices(M.shape[0],1)]
Xst=np.concatenate([iu(a[1]) for a in sites]); Xse=np.concatenate([iu(a[2]) for a in sites])
Xst=(Xst-Xst.mean())/Xst.std(); Xse=(Xse-Xse.mean())/Xse.std()
X=sm.add_constant(np.c_[Xst,Xse])
def fit_beta(y):
    yy=(y*(len(y)-1)+0.5)/len(y)  # squeeze 0/1 into (0,1)
    try:
        m=sm.GLM(yy,X,family=sm.families.Binomial(link=sm.families.links.Logit())).fit(scale='X2')  # quasi-binomial logit ~ beta-like
        return m.params[1],m.params[2],1-m.deviance/m.null_deviance
    except Exception: return np.nan,np.nan,np.nan
y0=np.concatenate([iu(a[0]) for a in sites])
b_st,b_se,pr2=fit_beta(y0)
print(f"n pairs={len(y0)} | Bray-Curtis=1.0 비율 {100*(y0>=0.999).mean():.0f}%")
print("\n=== 포화모델 (logit-link GLM, 구조+분광 동시) ===")
print(f"  구조 거리 계수 = {b_st:+.3f} | 분광 거리 계수 = {b_se:+.3f}  (둘 다 양수 = 거리↑ -> 비유사도↑)")
print(f"  pseudo-R2 = {pr2:.3f}")
# within-site permutation p for each coefficient
rng=np.random.RandomState(11); B=999; g_st=g_se=0
for _ in range(B):
    yp=np.concatenate([iu(a[0][np.ix_(p:=rng.permutation(a[3]),p)]) for a in sites])
    bs,be,_=fit_beta(yp)
    if bs>=b_st: g_st+=1
    if be>=b_se: g_se+=1
print(f"\n  within-site 순열 p: 구조 p={(g_st+1)/(B+1):.3f} | 분광 p={(g_se+1)/(B+1):.3f}")
print("  => 포화를 logit로 흡수하고 구조·분광을 동시통제해도 둘 다 유의하면 강건")
