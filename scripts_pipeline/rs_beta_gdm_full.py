# #4 GDM-style: monotone I-spline-like (hinge, non-negative coefs) transforms of structural &
# spectral distance -> Bray-Curtis; deviance explained + relative importance; within-site perm.
import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.spatial.distance import pdist, squareform
from sklearn.linear_model import LinearRegression
VS="E:/neon_lidar/vegetation_structure"; D="NEON_v2/data/"
mt=pd.read_csv(f"{VS}/vst_mappingandtagging.csv",usecols=['individualID','plotID','taxonID','siteID'],low_memory=False)
ai=pd.read_csv(f"{VS}/vst_apparentindividual.csv",usecols=['individualID','plantStatus','stemDiameter'],low_memory=False)
df=ai.merge(mt,on='individualID',how='left'); df=df[df.plantStatus.astype(str).str.contains('Live',na=False)]
df=df[df.stemDiameter>=10]; df=df[df.taxonID.notna() & ~df.taxonID.astype(str).str.contains('2PLANT|UNK',na=False)].drop_duplicates('individualID')
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
dst=np.concatenate([iu(a[1]) for a in sites]); dse=np.concatenate([iu(a[2]) for a in sites])
def ispline(x,k):  # monotone basis: x and hinges at quantile knots
    knots=np.quantile(x,np.linspace(0,1,k+2))[1:-1]
    return np.column_stack([x]+[np.maximum(0,x-kk) for kk in knots])
Bst=ispline(dst,3); Bse=ispline(dse,3); X=np.column_stack([Bst,Bse])
def fit(y):
    m=LinearRegression(positive=True).fit(X,y)  # non-neg coefs -> monotone increasing
    yh=m.predict(X); r2=1-((y-yh)**2).sum()/((y-y.mean())**2).sum()
    # relative importance = variance of each predictor's fitted transform
    fst=Bst@m.coef_[:Bst.shape[1]]; fse=Bse@m.coef_[Bst.shape[1]:]
    return r2, np.std(fst), np.std(fse)
y0=np.concatenate([iu(a[0]) for a in sites])
r2,ist,ise=fit(y0)
print(f"GDM-style (단조 spline): deviance/R² explained = {r2:.3f}")
print(f"  상대기여(적합 transform SD): 구조 {ist:.3f} | 분광 {ise:.3f} -> 구조 {100*ist/(ist+ise):.0f}% : 분광 {100*ise/(ist+ise):.0f}%")
rng=np.random.RandomState(11); B=499; ge=0
for _ in range(B):
    yp=np.concatenate([iu(a[0][np.ix_(p:=rng.permutation(a[3]),p)]) for a in sites])
    if fit(yp)[0]>=r2: ge+=1
print(f"  within-site 순열 p(모델 R²) = {(ge+1)/(B+1):.3f}")
