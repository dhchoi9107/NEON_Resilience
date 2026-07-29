# Within-site restricted-permutation Mantel: does structural/spectral dissimilarity track
# compositional dissimilarity (Bray-Curtis) plot-to-plot within sites?
import pandas as pd, numpy as np
from scipy.spatial.distance import pdist, squareform
VS="E:/neon_lidar/vegetation_structure"; D="NEON_v2/data/"
mt=pd.read_csv(f"{VS}/vst_mappingandtagging.csv",usecols=['individualID','plotID','taxonID','siteID'],low_memory=False)
ai=pd.read_csv(f"{VS}/vst_apparentindividual.csv",usecols=['individualID','plantStatus','stemDiameter'],low_memory=False)
df=ai.merge(mt,on='individualID',how='left')
df=df[df.plantStatus.astype(str).str.contains('Live',na=False)]; df=df[df.stemDiameter>=10]
df=df[df.taxonID.notna() & ~df.taxonID.astype(str).str.contains('2PLANT|UNK',na=False)].drop_duplicates('individualID')
f=pd.read_csv(D+"FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()
Scols=['VCI_mean','LAI_mean','Rugosity_mean','Vert_CV_mean','Canopy_Ht_mean','FHD_mean','Deep_Gap_mean']
Scols=[c for c in Scols if c in f.columns]; Pcols=[c for c in ['NDVI_mean','EVI_mean','ARVI_mean','SAVI_mean'] if c in f.columns]
for c in Scols+Pcols: f['z_'+c]=(f[c]-f[c].mean())/f[c].std()
# per-site distance matrices (species Bray-Curtis, structural & spectral Euclidean) on common plots
sites=[]
for s,g in f.groupby('siteID'):
    plots=g.plotID.tolist()
    sad=df[(df.siteID==s)&(df.plotID.isin(plots))].groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    common=[p for p in plots if p in sad.index]
    gg=g.set_index('plotID').loc[common]
    if len(common)<6 or gg[['z_'+c for c in Scols]].isna().any().any() or gg[['z_'+c for c in Pcols]].isna().any().any(): continue
    Dsp=squareform(pdist(sad.loc[common].values,metric='braycurtis'))
    Dst=squareform(pdist(gg[['z_'+c for c in Scols]].values))
    Dse=squareform(pdist(gg[['z_'+c for c in Pcols]].values))
    sites.append((Dsp,Dst,Dse,len(common)))
print(f"Mantel 사용 사이트: {len(sites)}")
iu=lambda M: M[np.triu_indices(M.shape[0],1)]
def pooled_r(perm_species=False,rng=None):
    xs_st=[];xs_se=[];ys=[]
    for Dsp,Dst,Dse,n in sites:
        if perm_species:
            p=rng.permutation(n); Dsp=Dsp[np.ix_(p,p)]
        ys.append(iu(Dsp)); xs_st.append(iu(Dst)); xs_se.append(iu(Dse))
    y=np.concatenate(ys); return np.corrcoef(np.concatenate(xs_st),y)[0,1], np.corrcoef(np.concatenate(xs_se),y)[0,1]
r_st,r_se=pooled_r()
rng=np.random.RandomState(11); B=999; ge_st=ge_se=0
for _ in range(B):
    a,b=pooled_r(True,rng)
    if a>=r_st: ge_st+=1
    if b>=r_se: ge_se+=1
print("\n=== Within-site restricted-permutation Mantel (999 perms) ===")
print(f"  구조 거리  ~ 조성 Bray-Curtis: Mantel r = {r_st:+.3f}, p = {(ge_st+1)/(B+1):.3f}")
print(f"  분광 거리  ~ 조성 Bray-Curtis: Mantel r = {r_se:+.3f}, p = {(ge_se+1)/(B+1):.3f}")
print("  (플롯을 사이트 내에서만 치환 -> 생물지리·사이트 구조 통제한 유효 p값)")
