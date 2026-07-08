"""
NEON_v2 STEP 5 — Temporal figures.
Builds per_year_v2 = LiDAR (all years) + VI (8 NEON-BRDF years) merged on plot-year.
  F07_timeseries_by_domain.png : domain-mean trajectories (VI shows 8 BRDF years only)
  F08_trend_distributions.png  : per-plot trend-slope distributions (net change direction)
  F09_coupling_heatmaps.png    : VI<->structure coupling across mean / SD / trend
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience"
V2=os.path.join(ROOT,"NEON_v2"); D,F,R=os.path.join(V2,"data"),os.path.join(V2,"figures"),os.path.join(V2,"results")
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=118,bbox_inches='tight'); plt.close(fig); print(" saved",n)

# ---- build per_year_v2 (LiDAR all yr + VI 8 BRDF yr) ----
LID=['Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV','Gini','VCI','FHD','LAI','Q95','Ht_Ratio']
VI=['NDVI','EVI','ARVI','SAVI']
lid=pd.read_csv(os.path.join(ROOT,"_ARCHIVE_v1","NEON_FINAL","data","per_year_RS_long.csv"))
lid=lid[['plotID','siteID','year']+LID]
vi=pd.read_csv(os.path.join(D,"plot_vi_neon_brdf.csv"))
vi=vi[['plotID','year']+[f"{v}_mean" for v in VI]].rename(columns={f"{v}_mean":v for v in VI})
py=lid.merge(vi,on=['plotID','year'],how='left')
DOM={"HARV":"D01","BART":"D01","SCBI":"D02","SERC":"D02","BLAN":"D02","GRSM":"D07","MLBS":"D07","ORNL":"D07",
 "JERC":"D03","OSBS":"D03","TALL":"D08","UNDE":"D05","STEI":"D05","TREE":"D05","WREF":"D16","ABBY":"D16",
 "SOAP":"D17","TEAK":"D17","RMNP":"D10"}
py['domain']=py['siteID'].map(DOM)
py.to_csv(os.path.join(D,"per_year_v2.csv"),index=False)
print("per_year_v2:",py.shape,"| VI rows(non-NaN NDVI):",py['NDVI'].notna().sum())

IDX=[('NDVI','VI'),('EVI','VI'),('SAVI','VI'),('ARVI','VI'),('LAI','LiDAR'),('FHD','LiDAR'),('Canopy_Ht','LiDAR'),('Deep_Gap','LiDAR')]

# ---- F07: domain-mean trajectories ----
doms=sorted(py['domain'].dropna().unique()); cmap=plt.cm.tab10(np.linspace(0,1,len(doms)))
fig,axes=plt.subplots(2,4,figsize=(22,10))
for ax,(idx,kind) in zip(axes.ravel(),IDX):
    for dm,c in zip(doms,cmap):
        d=py[py['domain']==dm][['year',idx]].dropna()
        if d.empty: continue
        cnt=d.groupby('year').size(); gy=d.groupby('year')[idx].mean(); gy=gy[cnt>=3]
        if len(gy)>=2: ax.plot(gy.index,gy.values,'-o',color=c,ms=3,lw=1,alpha=.8,label=dm)
    dd=py[['year',idx]].dropna()
    if len(dd)>5:
        sl,ic,r,p,se=st.linregress(dd['year'],dd[idx]); xx=np.array([dd['year'].min(),dd['year'].max()])
        ax.plot(xx,ic+sl*xx,'k--',lw=2.5,label=f"overall {sl:+.4f}/yr p={p:.0e}")
    tag=' (8 BRDF yr)' if kind=='VI' else ''
    ax.set_title(f"{idx} ({kind}){tag}",fontsize=12);ax.set_xlabel("year");ax.grid(alpha=.3)
axes.ravel()[0].legend(fontsize=7,ncol=2)
fig.suptitle("F07. RS time series by domain (VI = NEON .002 BRDF, 8 years; LiDAR all years; dashed=overall trend)",fontsize=14)
fig.tight_layout(); save(fig,"F07_timeseries_by_domain.png")

# ---- F08: per-plot trend slope distributions ----
fig,axes=plt.subplots(2,4,figsize=(20,9))
for ax,(idx,kind) in zip(axes.ravel(),IDX):
    slopes=[]
    for pid,g in py.groupby('plotID'):
        d=g[['year',idx]].dropna()
        if d['year'].nunique()>=3:
            sl=np.polyfit(d['year'],d[idx],1)[0]; slopes.append(sl)
    slopes=np.array(slopes)
    if len(slopes):
        ax.hist(slopes,bins=30,color='#2e7d32' if kind=='VI' else '#c62828',alpha=.7)
        ax.axvline(0,color='k',ls='--');ax.axvline(np.median(slopes),color='blue',lw=2,label=f"median={np.median(slopes):+.4f}")
        frac_pos=(slopes>0).mean()
        ax.set_title(f"{idx} ({kind}) n={len(slopes)}  {frac_pos*100:.0f}% increasing",fontsize=11);ax.legend(fontsize=8)
    ax.set_xlabel("per-plot slope /yr");ax.grid(alpha=.3)
fig.suptitle("F08. Per-plot temporal trend distributions (>=3 years)",fontsize=14)
fig.tight_layout(); save(fig,"F08_trend_distributions.png")

# ---- F09: VI<->structure coupling (pooled features) ----
pooled=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
pooled=pooled[pooled['sample_coverage']>=0.9]
STRUCT=['LAI','FHD','VCI','Canopy_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Gini']
fig,axes=plt.subplots(1,3,figsize=(20,6))
for ax,feat in zip(axes,['mean','sd','trend']):
    M=np.full((len(VI),len(STRUCT)),np.nan)
    for i,v in enumerate(VI):
        for j,s in enumerate(STRUCT):
            a=pooled.get(f"{v}_{feat}"); b=pooled.get(f"{s}_{feat}")
            if a is None or b is None: continue
            d=pd.concat([a,b],axis=1).dropna()
            if len(d)>20: M[i,j]=d.iloc[:,0].corr(d.iloc[:,1])
    im=ax.imshow(M,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
    ax.set_xticks(range(len(STRUCT)));ax.set_xticklabels(STRUCT,rotation=45,ha='right',fontsize=9)
    ax.set_yticks(range(len(VI)));ax.set_yticklabels(VI,fontsize=10)
    ax.set_title(f"VI–structure r ({feat})")
    for i in range(len(VI)):
        for j in range(len(STRUCT)):
            if np.isfinite(M[i,j]): ax.text(j,i,f"{M[i,j]:.2f}",ha='center',va='center',fontsize=7)
    fig.colorbar(im,ax=ax,fraction=.046)
fig.suptitle("F09. VI–structure coupling across temporal features (mean / interannual SD / trend)",fontsize=14)
fig.tight_layout(); save(fig,"F09_coupling_heatmaps.png")
print("DONE temporal figures")
