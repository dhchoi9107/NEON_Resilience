"""
Obj 2 + Obj 3 figures (from FINAL_v2_obj23.csv).
  G01_dhi_forest.png       : diversity ~ Sentinel DHI components
  G02_dhi_scatter.png      : Hill q1/q2 vs DHI cumulative/min/variation
  G03_disturbance_map.png  : plot locations colored by disturbance type
  G04_coupling_disturbed.png: RS<->diversity coupling, disturbed vs undisturbed
  G05_diversity_by_dist.png: diversity vs years-since-disturbance & by type
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=118,bbox_inches='tight'); plt.close(fig); print(" saved",n)

m=pd.read_csv(os.path.join(D,"FINAL_v2_obj23.csv"))
mc=m[m['sample_coverage']>=0.9].copy()

# G01: DHI forest
o3=pd.read_csv(os.path.join(R,"obj23_dhi_coeff.csv"))
if len(o3):
    resps=o3['response'].unique()
    fig,ax=plt.subplots(figsize=(9,6))
    yl=[]; i=0
    for resp in resps:
        for _,x in o3[o3.response==resp].iterrows():
            ax.errorbar(x['beta'],i,xerr=1.96*x['se'],fmt='o',
                        color='#1565c0' if x['p']<0.05 else '#90caf9',ms=7)
            yl.append(f"{resp} ~ {x['dhi'].replace('dhi_','').replace('_mean','').replace('_',' ')}"); i+=1
    ax.set_yticks(range(len(yl)));ax.set_yticklabels(yl,fontsize=8)
    ax.axvline(0,color='k',ls='--');ax.set_xlabel("beta (z)")
    ax.set_title("G01. Diversity ~ Sentinel-2 DHI (Obj 3 productivity)")
    fig.tight_layout(); save(fig,"G01_dhi_forest.png")

# G02: DHI scatter
comps=[('dhi_cum_mean','DHI cumulative'),('dhi_min_mean','DHI minimum'),('dhi_var_mean','DHI variation')]
fig,axes=plt.subplots(2,3,figsize=(15,8.5))
for j,(resp,rl) in enumerate([('Hill_q1','Hill q1'),('Hill_q2','Hill q2')]):
    for k,(c,cl) in enumerate(comps):
        ax=axes[j,k]; d=mc[[c,resp]].dropna()
        if len(d)<10: ax.axis('off'); continue
        ax.scatter(d[c],d[resp],s=10,alpha=.3,color='#1565c0')
        sl,ic,r,p,se=st.linregress(d[c],d[resp]);xx=np.linspace(d[c].min(),d[c].max(),30)
        ax.plot(xx,ic+sl*xx,'k-',lw=1.8)
        star='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        ax.set_title(f"{rl} vs {cl}  r={r:+.2f}{star}",fontsize=10);ax.set_xlabel(cl);ax.set_ylabel(rl)
fig.suptitle("G02. Pooled diversity vs Sentinel-2 DHI components",fontsize=14)
fig.tight_layout(); save(fig,"G02_dhi_scatter.png")

# G03: disturbance map
fig,ax=plt.subplots(figsize=(13,7))
col={'none':'#cfd8dc','fire':'#d32f2f','loss':'#f57c00','spectral':'#7b1fa2'}
# approximate lon/lat from siteID centroid not available -> use any easting/northing? use domain jitter
from pyproj import Transformer
ppy=pd.read_csv('E:/neon_lidar/vegetation_structure/vst_perplotperyear.csv',low_memory=False)
co=ppy.dropna(subset=['easting','northing','utmZone']).drop_duplicates('plotID')
def ll(r):
    z=int(str(r.utmZone)[:2]); tr=Transformer.from_crs(f"EPSG:326{z:02d}","EPSG:4326",always_xy=True)
    return tr.transform(r.easting,r.northing)
co[['lon','lat']]=co.apply(lambda r: pd.Series(ll(r)),axis=1)
mm=m.merge(co[['plotID','lon','lat']],on='plotID',how='left')
mm=mm[(mm.lon.between(-128,-65))&(mm.lat.between(24,50))]
for t,c in col.items():
    s=mm[mm.dist_type==t]
    ax.scatter(s.lon,s.lat,s=18,c=c,label=f"{t} ({len(s)})",alpha=.7,edgecolors='none')
ax.set_xlabel("lon");ax.set_ylabel("lat");ax.legend(title="disturbance");ax.grid(alpha=.3)
ax.set_title("G03. NEON plots by disturbance source (MTBS fire / Hansen loss / S2 spectral)")
fig.tight_layout(); save(fig,"G03_disturbance_map.png")

# G04: coupling disturbed vs undisturbed
o2=pd.read_csv(os.path.join(R,"obj23_disturbance.csv")) if os.path.exists(os.path.join(R,"obj23_disturbance.csv")) else pd.DataFrame()
pairs=[('SAVI_mean','Hill_q1'),('Deep_Gap_mean','Hill_q1'),('FHD_mean','Hill_q1'),('dhi_cum_mean','Hill_q1')]
fig,axes=plt.subplots(1,4,figsize=(20,5))
for ax,(rs,resp) in zip(axes,pairs):
    for dz,c,lab in [(0,'#2e7d32','undisturbed'),(1,'#d32f2f','disturbed')]:
        d=mc[mc.disturbed==dz][[rs,resp]].dropna()
        if len(d)<8: continue
        ax.scatter(d[rs],d[resp],s=12,alpha=.3,color=c,label=f"{lab} (n={len(d)})")
        sl,ic,r,p,se=st.linregress(d[rs],d[resp]);xx=np.linspace(d[rs].min(),d[rs].max(),30)
        ax.plot(xx,ic+sl*xx,'-',color=c,lw=2)
    ax.set_xlabel(rs);ax.set_ylabel(resp);ax.legend(fontsize=8)
    ax.set_title(f"{resp} ~ {rs}")
fig.suptitle("G04. RS<->taxonomic coupling: disturbed vs undisturbed (Obj 2 moderation)",fontsize=14)
fig.tight_layout(); save(fig,"G04_coupling_disturbed.png")

# G05: diversity by years-since-disturbance and by type
fig,axes=plt.subplots(1,2,figsize=(15,5.5))
d=mc[mc.disturbed==1][['years_since_dist','Hill_q1']].dropna()
axes[0].scatter(d['years_since_dist'],d['Hill_q1'],s=14,alpha=.4,color='#d32f2f')
if len(d)>5:
    sl,ic,r,p,se=st.linregress(d['years_since_dist'],d['Hill_q1']);xx=np.linspace(0,d['years_since_dist'].max(),20)
    axes[0].plot(xx,ic+sl*xx,'k-',lw=2,label=f"r={r:+.2f} p={p:.2g}")
axes[0].set_xlabel("years since disturbance");axes[0].set_ylabel("Hill q1");axes[0].legend()
axes[0].set_title("Recovery: diversity vs time since disturbance")
order=['none','spectral','loss','fire']
data=[mc[mc.dist_type==t]['Hill_q1'].dropna() for t in order]
axes[1].boxplot([x for x in data if len(x)],labels=[t for t,x in zip(order,data) if len(x)],showfliers=False)
axes[1].set_ylabel("Hill q1");axes[1].set_title("Diversity by disturbance type")
fig.suptitle("G05. Diversity, disturbance recency & type (Obj 2)",fontsize=14)
fig.tight_layout(); save(fig,"G05_diversity_by_dist.png")
print("DONE obj23 figures")
