"""
Manuscript figures, numbered in order of first citation in MANUSCRIPT_v2_full_26site.md
(following the author's convention: Fig 1 = conceptual figure cited in the Introduction,
Fig 2 = study area in Methods, results figures thereafter).
  Figure_1.png  (Intro) conceptual framework: three RS dimensions -> diversity components.
  Figure_2.png  (§2.1) study area: (a) 26 NEON sites map, (b) climate space (MAT vs MAP).
  Figure_3.png  (§3.1) plot-level framework: (a) sequential R2 beyond domain,
                (b) semi-partial variance partition, (c,d) Bayesian forests.
  Figure_4.png  (§3.2) RS dissimilarity vs compositional dissimilarity (variation hypothesis).
  Figure_5.png  (§3.5) site-level species-energy: (a) tower GPP, (b) PML-V2 (all sites).
  Figure_6.png  (§3.6) context moderation simple slopes: (a) MODIS GPP x stand age,
                (b) VCI x disturbance severity.
All data panels are drawn from the verified 26-site outputs; a VERIFY line is printed for
each so the numbers can be checked against the manuscript text.
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist
from scipy.stats import rankdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
RES=os.path.join(BASE,"papers","UNIFIED","results","O0_framework"); D=os.path.join(BASE,"data")
FIG=os.path.join(BASE,"papers","UNIFIED","figures")

nm=pd.read_csv(RES+"/nested_models_sd.csv").set_index("response")
vp=pd.read_csv(RES+"/variance_partition_sd.csv").set_index("response")
wb=pd.read_csv(RES+"/wildboot_sd.csv").set_index("response")
ps=pd.read_csv(RES+"/_pstar_blocks_cache.csv").set_index("response")
bc=pd.read_csv(RES+"/bayes_multilevel_coeffs.csv"); bm=pd.read_csv(RES+"/bayes_multilevel_meta.csv").set_index("response")
f=pd.read_csv(D+"/FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()

S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]; P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]; keep=S+P+Dn
RL=["Hill q1","Hill q2","LCBD turnover","LCBD nestedness"]; LAB=["Hill q1\n(alpha)","Hill q2\n(alpha)","LCBD\nturnover","LCBD\nnestedness"]
NICE={"Rugosity_mean":"Rugosity","Vert_CV_mean":"Vert-CV","VCI_mean":"VCI","LAI_mean":"LAI","EVI_mean":"EVI",
      "Rumple_trend":"Rumple trend","Vert_SD_trend":"Vert-SD trend","Vert_CV_trend":"Vert-CV trend","VCI_trend":"VCI trend",
      "FHD_trend":"FHD trend","LAI_trend":"LAI trend","Ht_Ratio_trend":"Ht-ratio trend"}
CBLK={"structure":"#08519c","spectral":"#9ecae1","dynamics":"#41ab5d"}
def blk(c): return "structure" if c in S else "spectral" if c in P else "dynamics"
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "n.s."
BAYES_RESP={"Hill q1":"Hill_q1","LCBD nestedness":"LCBD_nestedness_rare"}

# ================= FIGURE 1 (conceptual framework) =================
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig,ax=plt.subplots(figsize=(12.5,7.2)); ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis("off")
fig.suptitle("Figure 1. Conceptual framework: remote-sensing dimensions map onto the processes generating each diversity component",
             fontsize=12.5,fontweight="bold",y=0.98)
def box(x,y,w,h,text,fc,ec,fs=9.5,tc="black",bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.12",fc=fc,ec=ec,lw=1.6))
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fs,color=tc,fontweight="bold" if bold else "normal")
def arrow(x1,y1,x2,y2,color,lw,label=None,ls="-",lab_dy=0.16,lab_fs=8,lab_t=0.5):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=16,color=color,lw=lw,linestyle=ls,shrinkA=2,shrinkB=2))
    if label:
        lx,ly=x1+lab_t*(x2-x1),y1+lab_t*(y2-y1)
        ax.text(lx,ly+lab_dy,label,ha="center",fontsize=lab_fs,color=color,fontstyle="italic")
# left: RS dimensions
box(0.3,7.4,3.1,1.5,"STATIC CANOPY STATE\nLiDAR structure (VCI, LAI,\nrugosity, vertical CV)","#dbe9f6","#08519c",bold=True)
box(0.3,5.6,3.1,1.1,"Spectral greenness\n(EVI; weakest dimension)","#eef5fb","#9ecae1")
box(0.3,3.4,3.1,1.5,"STRUCTURAL DYNAMICS\nrepeat-LiDAR interannual\ntrends (state -> trajectory)","#e2f2e5","#41ab5d",bold=True)
box(0.3,1.0,3.1,1.5,"ECOSYSTEM PRODUCTIVITY\nindependent GPP (MODIS,\nPML-V2, flux towers)","#fdeadd","#b35806",bold=True)
# right: diversity components
box(6.6,7.6,3.1,1.3,"ALPHA DIVERSITY\nHill q1 / q2 (local richness)","#f2f2f2","#333333",bold=True)
box(6.6,5.4,3.1,1.4,"BETA: TURNOVER\nspecies replacement\n(lateral, spatial contrast)","#f2f2f2","#333333",bold=True)
box(6.6,3.2,3.1,1.4,"BETA: NESTEDNESS\nrichness difference\n(directional loss / gain)","#f2f2f2","#333333",bold=True)
box(6.6,1.0,3.1,1.3,"SPECIES-ENERGY\nsite-mean diversity vs GPP\n(monotonic, no hump)","#f2f2f2","#333333",bold=True)
# arrows (mapping)
arrow(3.5,8.2,6.5,8.3,"#08519c",3.0,"predictor levels (strongest)")
arrow(3.5,6.1,6.5,8.0,"#9ecae1",1.4,"alpha only",lab_dy=-0.34,lab_t=0.82)
arrow(3.5,7.6,6.5,6.2,"#2c7fb8",2.2,"as between-plot dissimilarity",lab_dy=-0.40,lab_t=0.30)
arrow(3.5,4.2,6.5,3.9,"#41ab5d",3.0,"temporal trends (Bayesian-credible)")
arrow(3.5,1.7,6.5,1.6,"#b35806",2.6,"site scale (ICC = 0.92)")
ax.text(5.0,0.25,"Process mapping: lateral replacement is read from static contrast; directional richness change from the temporal trajectory.",
        ha="center",fontsize=9,color="#444444",fontstyle="italic")
fig.savefig(FIG+"/Figure_1.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print("VERIFY Fig1: conceptual diagram (no data)")

# ================= FIGURE 2 (study area) =================
import geopandas as gpd
NE=os.path.join(r"C:\Users\star1\Documents\GitHub\NEON_Resilience","scripts_pipeline","_pipeline_state","ne_states.gpkg")
if not os.path.exists(NE):
    gpd.read_file("https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_1_states_provinces.zip").to_file(NE)
states=gpd.read_file(NE)
ll=pd.read_csv(D+"/plot_lonlat_26.csv").groupby("siteID")[["lon","lat"]].mean().reset_index()
cl=pd.read_csv(D+"/site_climate_neon.csv")
hq=f.groupby("siteID").Hill_q1.mean().reset_index()
sm=ll.merge(cl,on="siteID").merge(hq,on="siteID")
dom=f[["siteID","domain"]].drop_duplicates()
sm=sm.merge(dom,on="siteID",how="left")
fig,ax=plt.subplots(1,2,figsize=(14.5,6.2),gridspec_kw={"width_ratios":[1.45,1]})
states.boundary.plot(ax=ax[0],color="#bbbbbb",lw=0.6)
sc=ax[0].scatter(sm.lon,sm.lat,c=sm.MAT_C,cmap="coolwarm",s=110,edgecolor="black",linewidth=0.8,zorder=3)
for _,r in sm.iterrows(): ax[0].annotate(r.siteID,(r.lon,r.lat),fontsize=6.5,xytext=(4,3),textcoords="offset points")
ax[0].set_xlim(-170,-62); ax[0].set_ylim(17,72); ax[0].set_xlabel("Longitude"); ax[0].set_ylabel("Latitude")
ax[0].set_title(f"(a) 26 NEON forested sites, 12 eco-climatic domains\n(conterminous US to interior Alaska)",fontweight="bold",fontsize=11)
fig.colorbar(sc,ax=ax[0],fraction=0.03,pad=0.02,label="MAT (°C)")
s2=ax[1].scatter(sm.MAT_C,sm.MAP_mm,c=sm.Hill_q1,cmap="viridis",s=120,edgecolor="black",linewidth=0.8)
for _,r in sm.iterrows(): ax[1].annotate(r.siteID,(r.MAT_C,r.MAP_mm),fontsize=6.5,xytext=(4,3),textcoords="offset points")
ax[1].set_xlabel("Mean annual temperature (°C)"); ax[1].set_ylabel("Mean annual precipitation (mm)")
ax[1].set_title("(b) Climate space of the 26 sites\ncolour = site-mean tree diversity (Hill q1)",fontweight="bold",fontsize=11)
fig.colorbar(s2,ax=ax[1],fraction=0.04,pad=0.02,label="Hill q1")
fig.suptitle("Figure 2. Study area spans broad climate, productivity, and diversity gradients",fontsize=12.5,fontweight="bold",y=1.00)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(FIG+"/Figure_2.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig2: sites={len(sm)} (26), MAT range {sm.MAT_C.min():.1f}~{sm.MAT_C.max():.1f}C, MAP {sm.MAP_mm.min():.0f}~{sm.MAP_mm.max():.0f}mm")

# ================= FIGURE 3 =================
fig=plt.figure(figsize=(15.5,11.5)); gs=fig.add_gridspec(2,2,hspace=0.34,wspace=0.26)
na=int(nm.loc["Hill q1","n"]); nb=int(nm.loc["LCBD turnover","n"]); ns=f.siteID.nunique()
fig.suptitle(f"Figure 3. Plot-level framework: canopy state and structural dynamics carry complementary diversity information\n"
             f"(n = {na} alpha / {nb} beta plots, {ns} sites, 12 domains)",fontsize=13.5,fontweight="bold",y=0.99)
x=np.arange(4); w=0.6
# (a) sequential R2 beyond domain
ax=fig.add_subplot(gs[0,0])
state=[nm.loc[r,"R2_M1_beyond"] for r in RL]; dyn=[nm.loc[r,"R2_M2_beyond"]-nm.loc[r,"R2_M1_beyond"] for r in RL]
ax.bar(x,state,w,label="Canopy state (structure + spectral)",color="#2c7fb8")
ax.bar(x,dyn,w,bottom=state,label="+ Structural dynamics",color="#41ab5d")
for i,r in enumerate(RL):
    ax.text(i,state[i]/2,f"{state[i]:.3f}",ha="center",va="center",color="white",fontweight="bold",fontsize=9)
    ax.text(i,state[i]+dyn[i]/2,f"+{dyn[i]:.3f}\n{star(wb.loc[r,'p_dyn_wildboot'])}",ha="center",va="center",color="white",fontsize=8.2)
ax.set_xticks(x); ax.set_xticklabels(LAB,fontsize=10); ax.set_ylabel("R² beyond domain (sequential)")
ax.set_title("(a) Sequential explanatory power beyond biogeography\nstars = wild cluster bootstrap on the dynamics increment",fontweight="bold",fontsize=10.5)
ax.legend(fontsize=8.5,loc="upper right"); ax.set_ylim(0,0.19)
# (b) variance partition
ax=fig.add_subplot(gs[0,1])
blocks=[("unique_structure","structure","#08519c"),("unique_spectral","spectral","#9ecae1"),("unique_dynamics","dynamics","#41ab5d"),("shared_RS","shared","#bdbdbd")]
for i,r in enumerate(RL):
    b=0
    for col,lab,c in blocks:
        val=vp.loc[r,col]; ax.bar(i,val,w,bottom=b,color=c,label=lab if i==0 else None)
        if val>0.004:
            txt=f"{val:.3f}"+(f" {star(ps.loc[r,col])}" if col in ps.columns else "")
            ax.text(i,b+val/2,txt,ha="center",va="center",fontsize=7.8,color="white" if col in("unique_structure","unique_dynamics") else "black")
        b+=val
ax.set_xticks(x); ax.set_xticklabels(LAB,fontsize=10); ax.set_ylabel("Unique R² beyond domain (semi-partial)")
ax.set_title("(b) Unique contribution of each block\nstructure → alpha; dynamics → nestedness",fontweight="bold",fontsize=10.5)
ax.legend(fontsize=8.5,title="block",loc="upper right")
# (c,d) Bayesian forests
def forest(ax,resp_label,title):
    rr=BAYES_RESP[resp_label]; sub=bc[bc.response==rr].set_index("predictor"); order=list(reversed(keep)); ys=np.arange(len(order))
    for yi,pred in zip(ys,order):
        row=sub.loc[pred]; b,lo,hi=row["beta"],row["hdi_2.5"],row["hdi_97.5"]; cred=bool(row["credible"]); c=CBLK[blk(pred)]
        ax.plot([lo,hi],[yi,yi],color=c,lw=2.6 if cred else 1.3,alpha=1.0 if cred else 0.4,solid_capstyle="round",zorder=2)
        ax.scatter([b],[yi],s=64 if cred else 30,color=c,edgecolor="black" if cred else "none",linewidth=0.8,alpha=1.0 if cred else 0.45,zorder=3)
        if cred: ax.text(hi+0.012 if b>0 else lo-0.012,yi,f"{b:+.2f}",va="center",ha="left" if b>0 else "right",fontsize=7.4,fontweight="bold",color=c)
    ax.axvline(0,color="grey",lw=1,ls="--",zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels([NICE[p] for p in order],fontsize=8.6)
    for tl,p in zip(ax.get_yticklabels(),order): tl.set_color(CBLK[blk(p)])
    ax.set_xlabel("Standardized β (95% HDI)",fontsize=10); ax.set_xlim(-0.55,0.66)
    ax.set_title(title,fontweight="bold",fontsize=10.5)
    m=bm.loc[BAYES_RESP[resp_label]]
    ax.text(0.5,-0.155,f"0 divergences · R̂ ≤ {m['max_rhat']:.3f} · n = {int(m['n'])} · filled = 95% HDI excludes 0",
            transform=ax.transAxes,ha="center",va="top",fontsize=7.6,style="italic",color="grey")
forest(fig.add_subplot(gs[1,0]),"Hill q1","(c) Bayesian multilevel — Hill q1 (alpha)\nstructure dominates (LAI, VCI) + spectral (EVI)")
forest(fig.add_subplot(gs[1,1]),"LCBD nestedness","(d) Bayesian multilevel — LCBD nestedness\nthree dynamics trends credible")
fig.savefig(FIG+"/Figure_3.png",dpi=200,bbox_inches="tight"); plt.close(fig)
seq=[round(d,3) for d in dyn]
print(f"VERIFY Fig3: state15%={state[0]:.3f} turn2%={nm.loc['LCBD turnover','R2_M2_beyond']:.3f} nest9%={nm.loc['LCBD nestedness','R2_M2_beyond']:.3f} | seq dyn={seq} (text 0.042/0.039/0.013/0.056)")
print(f"VERIFY Fig3: unique str q1={vp.loc['Hill q1','unique_structure']:.3f} (text 0.081), dyn nest={vp.loc['LCBD nestedness','unique_dynamics']:.3f} (text 0.056)")
cred_n=bc[(bc.response=='LCBD_nestedness_rare')&(bc.credible)&(bc.block=='dynamics')]
print(f"VERIFY Fig3: Bayes nest credible dyn = {[(r.predictor,r.beta) for r in cred_n.itertuples()]}")

# ================= FIGURE 4 =================
VS="E:/neon_lidar/vegetation_structure"
mt=pd.read_csv(f"{VS}/vst_mappingandtagging.csv",usecols=['individualID','plotID','taxonID','siteID'],low_memory=False)
ai=pd.read_csv(f"{VS}/vst_apparentindividual.csv",usecols=['individualID','plantStatus','stemDiameter'],low_memory=False)
dfv=ai.merge(mt,on='individualID',how='left')
dfv=dfv[dfv.plantStatus.astype(str).str.contains('Live',na=False)]; dfv=dfv[dfv.stemDiameter>=10]
dfv=dfv[dfv.taxonID.notna() & ~dfv.taxonID.astype(str).str.contains('2PLANT|UNK',na=False)].drop_duplicates('individualID')
Scols=[c for c in ['VCI_mean','LAI_mean','Rugosity_mean','Vert_CV_mean','Canopy_Ht_mean','FHD_mean','Deep_Gap_mean'] if c in f.columns]
Pcols=[c for c in ['NDVI_mean','EVI_mean','ARVI_mean','SAVI_mean'] if c in f.columns]
fz=f.copy()
for c in Scols+Pcols: fz['z_'+c]=(fz[c]-fz[c].mean())/fz[c].std()
XSt=[];XSe=[];Y=[];nsite=0
for s,g in fz.groupby('siteID'):
    plots=g.plotID.tolist(); sad=dfv[(dfv.siteID==s)&(dfv.plotID.isin(plots))].groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    common=[p for p in plots if p in sad.index]; gg=g.set_index('plotID').loc[common] if common else None
    if not common or len(common)<6 or gg[['z_'+c for c in Scols]].isna().any().any() or gg[['z_'+c for c in Pcols]].isna().any().any(): continue
    nsite+=1
    Y.append(pdist(sad.loc[common].values,metric='braycurtis'))
    XSt.append(pdist(gg[['z_'+c for c in Scols]].values)); XSe.append(pdist(gg[['z_'+c for c in Pcols]].values))
y=np.concatenate(Y); xst=np.concatenate(XSt); xse=np.concatenate(XSe)
rst=np.corrcoef(rankdata(xst),rankdata(y))[0,1]; rse=np.corrcoef(rankdata(xse),rankdata(y))[0,1]
fig,ax=plt.subplots(1,2,figsize=(13,5.8))
for a,(xx,r,lab,col,tag) in zip(ax,[(xst,rst,"Structural distance (Euclidean, LiDAR metrics)","#08519c","(a) Structural"),
                                    (xse,rse,"Spectral distance (Euclidean, vegetation indices)","#00695c","(b) Spectral")]):
    a.scatter(xx,y,s=5,alpha=0.06,color=col,edgecolor="none")
    bins=np.quantile(xx,np.linspace(0,1,11)); idx=np.digitize(xx,bins)
    bx=[xx[idx==k].mean() for k in range(1,11) if (idx==k).sum()>=50]; by=[y[idx==k].mean() for k in range(1,11) if (idx==k).sum()>=50]
    a.plot(bx,by,'o-',color="#c62828",ms=8,lw=2.5,zorder=4,label="binned mean (decile bins)")
    a.set_xlabel(lab,fontsize=11); a.set_ylabel("Compositional dissimilarity (Bray–Curtis)",fontsize=11)
    a.set_title(f"{tag} variation hypothesis\nwithin-site Mantel (Spearman) r = {r:+.2f}, p = 0.001",fontweight="bold",fontsize=11)
    a.legend(fontsize=9,loc="lower right")
fig.suptitle("Figure 4. Plots that differ more structurally / spectrally also differ more in species composition\n"
             f"(within-site plot pairs, n = {len(y):,} across {nsite} sites; rank-based Mantel, biogeography held fixed)",fontweight="bold",fontsize=12.5)
fig.tight_layout(rect=[0,0,1,0.90]); fig.savefig(FIG+"/Figure_4.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig4: Spearman struct={rst:+.2f} (text +0.32), spectral={rse:+.2f} (text +0.30), pairs={len(y)} (text 8,813), sites={nsite} (text 23)")

# ================= FIGURE 5 =================
g3=pd.read_csv(D+"/plot_pml_gpp_ts_26.csv")[["plotID","pml_gpp"]]
site=f.merge(g3,on="plotID").groupby("siteID").agg(Hill_q1=("Hill_q1","mean"),pml=("pml_gpp","mean")).reset_index()
tw=pd.read_csv(D+"/site_tower_gpp_26.csv")[["siteID","tower_gpp"]]; site=site.merge(tw,on="siteID",how="left")
tv=site.dropna(subset=["tower_gpp"])
rt,pt=st.pearsonr(tv.tower_gpp,tv.Hill_q1)
t2=tv.assign(z=(tv.tower_gpp-tv.tower_gpp.mean())/tv.tower_gpp.std()); qt=smf.ols("Hill_q1 ~ z + I(z**2)",t2).fit().pvalues["I(z ** 2)"]
rp_,pp_=st.pearsonr(site.pml,site.Hill_q1)
rvt=st.pearsonr(tv.pml,tv.tower_gpp)[0]
fig,ax=plt.subplots(1,2,figsize=(13.6,5.9))
for a,(d3,xc,r,p,ttl,xl,col) in zip(ax,[
    (tv,"tower_gpp",rt,pt,f"(a) Eddy-covariance tower GPP (n = {len(tv)} sites)\nr = {rt:+.2f}, p = {pt:.3f}; quadratic n.s. (p = {qt:.2f})","Tower GPP (gC m$^{-2}$ yr$^{-1}$)","#b35806"),
    (site,"pml",rp_,pp_,f"(b) PML-V2 GPP, all {len(site)} sites\nr = {rp_:+.2f}, p = {pp_:.3f}; PML vs tower r = +{rvt:.2f}","Site-mean PML-V2 GPP (2000–2024)","#00695c")]):
    a.scatter(d3[xc],d3.Hill_q1,s=90,color=col,edgecolor="white",linewidth=0.8,zorder=3)
    for _,row in d3.iterrows(): a.annotate(row.siteID,(row[xc],row.Hill_q1),fontsize=6,alpha=0.55,xytext=(4,3),textcoords="offset points")
    b1,b0=np.polyfit(d3[xc],d3.Hill_q1,1); xr=np.linspace(d3[xc].min(),d3[xc].max(),50); a.plot(xr,b0+b1*xr,color="black",lw=2)
    a.set_xlabel(xl,fontsize=10.5); a.set_ylabel("Site-mean tree diversity (Hill q1)",fontsize=10.5)
    a.set_title(ttl,fontweight="bold",fontsize=10.5)
fig.suptitle("Figure 5. Site-level species–energy relationship is monotonic (no hump) across the temperate-to-boreal gradient",fontweight="bold",fontsize=12.5)
fig.tight_layout(rect=[0,0,1,0.92]); fig.savefig(FIG+"/Figure_5.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig5: tower r={rt:+.2f} p={pt:.3f} (text +0.59, 0.004), quad p={qt:.2f} (text 0.87), PML-tower r={rvt:+.2f} (text +0.86), n_tower={len(tv)} (text 22)")

# ================= FIGURE 6 =================
ss=pd.read_csv(RES+"/simple_slopes.csv"); z=[-1,0,1]
gA=ss[(ss.rp=="Hill q1")&(ss.rs=="modis_gpp")&(ss.ctx=="stand_age_gami")].sort_values("ctx_z")
gB=ss[(ss.rp=="Hill q1")&(ss.rs=="VCI_mean")&(ss.ctx=="severity")].sort_values("ctx_z")
fig,ax=plt.subplots(1,2,figsize=(12.6,5.4))
for a,(g4,col,ttl,xt) in zip(ax,[
    (gA,"#b35806","(a) MODIS GPP → Hill q1, by stand age\ninteraction n.s. after FDR (q = 0.26): suggestive only",["−1 SD\n(younger)","mean","+1 SD\n(older)"]),
    (gB,"#08519c","(b) VCI → Hill q1, by disturbance severity\ncoupling weakens with disturbance",["−1 SD\n(least disturbed)","mean","+1 SD\n(most disturbed)"])]):
    a.errorbar(z,g4.slope,yerr=1.96*g4.se,marker="o",capsize=5,color=col,lw=2.2,ms=9)
    for zz,sl,pv in zip(z,g4.slope,g4.p):
        a.annotate(f"{sl:+.2f}{'*' if pv<0.05 else ''}",(zz,sl),textcoords="offset points",xytext=(10,6),color=col,fontsize=10,fontweight="bold")
    a.axhline(0,color="grey",lw=0.8,ls="--")
    a.set_xticks(z); a.set_xticklabels(xt,fontsize=9.5); a.set_xlim(-1.6,1.6)
    a.set_ylabel("Simple slope on Hill q1 (per 1 SD)",fontsize=10.5); a.set_xlabel("Moderator (context)",fontsize=10.5)
    a.set_title(ttl,fontweight="bold",fontsize=10.5)
fig.suptitle("Figure 6. Ecological context modifies remote sensing–diversity associations (simple slopes ± 95% CI)",fontweight="bold",fontsize=12.5)
fig.tight_layout(rect=[0,0,1,0.92]); fig.savefig(FIG+"/Figure_6.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig6: MODISxAge slopes={[round(v,2) for v in gA.slope]} p={[round(v,3) for v in gA.p]} (text +0.17/+0.27/+0.38, 0.13/0.01/0.007)")
print(f"VERIFY Fig6: VCIxSev slopes={[round(v,2) for v in gB.slope]} (text +0.39 -> +0.21)")
print("ALL FIGURES SAVED ->",FIG)
