import pandas as pd, numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
VS="E:/neon_lidar/vegetation_structure"; D="NEON_v2/data/"; FIG="NEON_v2/papers/UNIFIED/figures/"
mt=pd.read_csv(f"{VS}/vst_mappingandtagging.csv",usecols=['individualID','plotID','taxonID','siteID'],low_memory=False)
ai=pd.read_csv(f"{VS}/vst_apparentindividual.csv",usecols=['individualID','plantStatus','stemDiameter'],low_memory=False)
df=ai.merge(mt,on='individualID',how='left')
df=df[df.plantStatus.astype(str).str.contains('Live',na=False)]; df=df[df.stemDiameter>=10]
df=df[df.taxonID.notna() & ~df.taxonID.astype(str).str.contains('2PLANT|UNK',na=False)].drop_duplicates('individualID')
f=pd.read_csv(D+"FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()
Scols=[c for c in ['VCI_mean','LAI_mean','Rugosity_mean','Vert_CV_mean','Canopy_Ht_mean','FHD_mean','Deep_Gap_mean'] if c in f.columns]
Pcols=[c for c in ['NDVI_mean','EVI_mean','ARVI_mean','SAVI_mean'] if c in f.columns]
for c in Scols+Pcols: f['z_'+c]=(f[c]-f[c].mean())/f[c].std()
XSt=[];XSe=[];Y=[]
for s,g in f.groupby('siteID'):
    plots=g.plotID.tolist(); sad=df[(df.siteID==s)&(df.plotID.isin(plots))].groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    common=[p for p in plots if p in sad.index]; gg=g.set_index('plotID').loc[common]
    if len(common)<6 or gg[['z_'+c for c in Scols]].isna().any().any() or gg[['z_'+c for c in Pcols]].isna().any().any(): continue
    Y.append(pdist(sad.loc[common].values,metric='braycurtis'))
    XSt.append(pdist(gg[['z_'+c for c in Scols]].values)); XSe.append(pdist(gg[['z_'+c for c in Pcols]].values))
y=np.concatenate(Y); xst=np.concatenate(XSt); xse=np.concatenate(XSe)
rst=np.corrcoef(rankdata(xst),rankdata(y))[0,1]; rse=np.corrcoef(rankdata(xse),rankdata(y))[0,1]  # Spearman
fig,ax=plt.subplots(1,2,figsize=(13,5.8))
for a,(xx,r,lab,col) in zip(ax,[(xst,rst,"Structural distance (Euclidean, LiDAR metrics)","#08519c"),(xse,rse,"Spectral distance (Euclidean, VIs)","#00695c")]):
    a.scatter(xx,y,s=5,alpha=0.06,color=col,edgecolor="none")
    # binned means for readability
    bins=np.quantile(xx,np.linspace(0,1,11)); idx=np.digitize(xx,bins)
    bx=[xx[idx==k].mean() for k in range(1,11) if (idx==k).sum()>=50]; by=[y[idx==k].mean() for k in range(1,11) if (idx==k).sum()>=50]
    a.plot(bx,by,'o-',color="#c62828",ms=8,lw=2.5,zorder=4,label="binned mean (decile bins)")
    a.set_xlabel(lab,fontsize=11); a.set_ylabel("Compositional dissimilarity (Bray–Curtis)",fontsize=11)
    a.set_title(f"{'(A) Structural' if col=='#08519c' else '(B) Spectral'} variation hypothesis\nwithin-site Mantel (Spearman) r = {r:+.2f}, p = 0.001",fontweight="bold",fontsize=11)
    a.legend(fontsize=9,loc="lower right")
fig.suptitle("Plots that differ more structurally / spectrally also differ more in species composition\n(within-site plot pairs; rank-based Mantel, biogeography controlled; 9% of pairs share no species)",fontweight="bold",fontsize=12.5)
fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig(FIG+"F3_rs_variation_beta.png",dpi=200,bbox_inches="tight"); print("saved F3, pairs n=",len(y))
