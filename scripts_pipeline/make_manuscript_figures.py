"""
Publication-quality manuscript figures (GEB double-column, 180 mm width, Arial, despined),
numbered in order of first citation in MANUSCRIPT_v2_full_26site.md.
  Figure_1  (Intro)  conceptual framework
  Figure_2  (S2.1)   study area: (a) Albers map of 26 sites, (b) climate space
  Figure_3  (S3.1)   plot-level framework: (a) sequential R2, (b) variance partition,
                     (c,d) Bayesian multilevel forests
  Figure_4  (S3.2)   RS dissimilarity vs compositional dissimilarity (hexbin + decile means)
  Figure_5  (S3.5)   site-level species-energy: (a) tower GPP, (b) PML-V2
  Figure_6  (S3.6)   context moderation simple slopes
In-figure titles/captions are omitted (they live in the manuscript); panels carry only
letters and minimal labels. A VERIFY line is printed per figure against the text numbers.
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist
from scipy.stats import rankdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------- global publication style ----------------
plt.rcParams.update({
    "font.family":"Arial","font.size":8,"axes.titlesize":8.5,"axes.labelsize":8.5,
    "xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5,
    "axes.spines.top":False,"axes.spines.right":False,
    "axes.linewidth":0.8,"xtick.direction":"out","ytick.direction":"out",
    "xtick.major.width":0.8,"ytick.major.width":0.8,"xtick.major.size":3,"ytick.major.size":3,
    "savefig.dpi":400,"figure.dpi":120,"pdf.fonttype":42,"ps.fonttype":42})
W=7.09   # 180 mm text width
def panel(ax,letter,dx=-0.10,dy=1.04):
    ax.text(dx,dy,letter,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom",ha="left")

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
RL=["Hill q1","Hill q2","LCBD turnover","LCBD nestedness"]; LAB=["Hill q1","Hill q2","Turnover","Nestedness"]
NICE={"Rugosity_mean":"Rugosity","Vert_CV_mean":"Vertical CV","VCI_mean":"VCI","LAI_mean":"LAI","EVI_mean":"EVI",
      "Rumple_trend":"Rumple trend","Vert_SD_trend":"Vertical-SD trend","Vert_CV_trend":"Vertical-CV trend","VCI_trend":"VCI trend",
      "FHD_trend":"FHD trend","LAI_trend":"LAI trend","Ht_Ratio_trend":"Height-ratio trend"}
C_STR="#3a5a40"; C_SPE="#a3b18a"; C_DYN="#588157"; C_SHA="#c9c9c9"; C_PRO="#a98467"; C_ACC="#c1121f"   # Forest-earth palette (user-selected)
CBLK={"structure":C_STR,"spectral":C_SPE,"dynamics":C_DYN}
def blk(c): return "structure" if c in S else "spectral" if c in P else "dynamics"
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "ns"
BAYES_RESP={"Hill q1":"Hill_q1","LCBD nestedness":"LCBD_nestedness_rare"}

# ================= FIGURE 1 : conceptual framework =================
fig,ax=plt.subplots(figsize=(W,4.1)); ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis("off")
from matplotlib.colors import to_rgba
def cbox(x,y,w,h,lines,accent=None,fs=7.6):
    if accent:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.10",fc=to_rgba(accent,0.10),ec=accent,lw=1.3))
    else:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.10",fc="#f7f7f7",ec="#666666",lw=1.0))
    t0=lines[0]; rest="\n".join(lines[1:])
    ax.text(x+w/2,y+h-0.34,t0,ha="center",va="center",fontsize=fs+0.6,fontweight="bold",color="#111111")
    if rest: ax.text(x+w/2,y+(h-0.62)/2,rest,ha="center",va="center",fontsize=fs,color="#444444")
def carrow(x1,y1,x2,y2,color,lw,label,lab_t=0.5,lab_dy=0.22,curve=0.0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=11,color=color,lw=lw,
                                 shrinkA=3,shrinkB=3,connectionstyle=f"arc3,rad={curve}"))
    lx,ly=x1+lab_t*(x2-x1),y1+lab_t*(y2-y1)
    ax.text(lx,ly+lab_dy,label,ha="center",fontsize=6.8,color=color,fontstyle="italic")
LH=1.62; ys=[8.20,6.15,3.90,1.65]
cbox(0.25,ys[0],3.15,LH,["Static canopy structure","airborne LiDAR","(VCI, LAI, rugosity, vertical CV)"],accent=C_STR)
cbox(0.25,ys[1],3.15,LH,["Spectral greenness","BRDF-corrected VIs (EVI)","weakest dimension"],accent=C_SPE)
cbox(0.25,ys[2],3.15,LH,["Structural dynamics","repeat-LiDAR interannual trends","state → trajectory"],accent=C_DYN)
cbox(0.25,ys[3],3.15,LH,["Ecosystem productivity","independent GPP","(PML-V2 satellite, flux towers)"],accent=C_PRO)
cbox(6.60,ys[0],3.15,LH,["Alpha diversity","Hill q1 / q2","local richness"])
cbox(6.60,ys[1],3.15,LH,["Beta: turnover","species replacement","lateral, spatial contrast"])
cbox(6.60,ys[2],3.15,LH,["Beta: nestedness","richness difference","directional loss / gain"])
cbox(6.60,ys[3],3.15,LH,["Species–energy","site-mean diversity vs GPP","monotonic, no hump"])
yc=[y+LH/2 for y in ys]
carrow(3.45,yc[0]+0.25,6.55,yc[0]+0.25,C_STR,2.6,"predictor levels (strongest)",0.5,0.24)
carrow(3.45,yc[1]+0.30,6.55,yc[0]-0.30,"#8a9a6b",1.1,"alpha only",0.80,-0.34)
carrow(3.45,yc[0]-0.42,6.55,yc[1]+0.30,"#6b705c",1.8,"as between-plot dissimilarity",0.85,-0.36)
carrow(3.45,yc[2],6.55,yc[2],C_DYN,2.6,"temporal trends (Bayesian-credible)",0.5,0.24)
carrow(3.45,yc[3],6.55,yc[3],C_PRO,2.2,"site scale (ICC = 0.92)",0.5,0.24)
ax.text(5.0,0.28,"Lateral replacement is read from static contrast; directional richness change from the temporal trajectory.",
        ha="center",fontsize=7,color="#666666",fontstyle="italic")
fig.savefig(FIG+"/Figure_1.png",bbox_inches="tight"); plt.close(fig)
print("VERIFY Fig1: conceptual diagram (no data)")

# ================= FIGURE 2 : study area (Albers) =================
import geopandas as gpd
NE=os.path.join(r"C:\Users\star1\Documents\GitHub\NEON_Resilience","scripts_pipeline","_pipeline_state","ne_states.gpkg")
if not os.path.exists(NE):
    gpd.read_file("https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_1_states_provinces.zip").to_file(NE)
from adjustText import adjust_text
from shapely.geometry import box as _box
DOM=os.path.join(r"C:\Users\star1\Documents\GitHub\NEON_Resilience","scripts_pipeline","_pipeline_state","neon_domains.gpkg")
if not os.path.exists(DOM):
    gpd.read_file("https://www.neonscience.org/sites/default/files/NEONDomains_0.zip").to_file(DOM)
dom=gpd.read_file(DOM).to_crs(4326)
conus_dom=gpd.clip(dom,_box(-125.5,24.0,-66.5,49.7)).to_crs("EPSG:5070")   # NEON eco-climatic domains, CONUS
ak_dom=gpd.clip(dom,_box(-169.5,52.0,-129.5,71.8)).to_crs("EPSG:3338")      # Alaska domains
ll=pd.read_csv(D+"/plot_lonlat_26.csv").groupby("siteID")[["lon","lat"]].mean().reset_index()
cl=pd.read_csv(D+"/site_climate_neon.csv"); hq=f.groupby("siteID").Hill_q1.mean().reset_index()
sm=ll.merge(cl,on="siteID").merge(hq,on="siteID")
gsm=gpd.GeoDataFrame(sm,geometry=gpd.points_from_xy(sm.lon,sm.lat),crs=4326)
p48=gsm[gsm.lat<50].to_crs("EPSG:5070"); pak=gsm[gsm.lat>=50].to_crs("EPSG:3338")
import matplotlib.patheffects as pe
HALO=[pe.withStroke(linewidth=1.8,foreground="white")]
def autolab(axh,xs,ys,names,fs=4.8,colors=None):
    if colors is None: colors=["#333333"]*len(names)
    txts=[axh.text(x_,y_,n_,fontsize=fs,color=c_,fontweight="bold",path_effects=HALO)
          for x_,y_,n_,c_ in zip(xs,ys,names,colors)]
    adjust_text(txts,x=list(xs),y=list(ys),ax=axh,expand=(1.25,1.6),force_text=(0.4,0.6),
                arrowprops=dict(arrowstyle="-",color="#999999",lw=0.4))
# site labels coloured by forest type (conifer stem fraction from vst; >=75% needleleaf, <=25% broadleaf)
FTYPE=pd.read_csv(os.path.join(r"C:\Users\star1\Documents\GitHub\NEON_Resilience","scripts_pipeline","_pipeline_state","site_foresttype.csv"))
FTC={"Needleleaf":"#1b6535","Mixed":"#5a5a5a","Broadleaf":"#a4551e"}
sm=sm.merge(FTYPE[["siteID","forest_type"]],on="siteID",how="left")
sm["labcol"]=sm.forest_type.map(FTC).fillna("#333333")
# --- continuous bivariate colour: MAT x MAP, bilinear blend of four corner colours ---
_c00,_c10=np.array(to_rgba("#e8e8e8")[:3]),np.array(to_rgba("#c0392b")[:3])   # dry: cool -> warm
_c01,_c11=np.array(to_rgba("#4d9ec4")[:3]),np.array(to_rgba("#3f2a3d")[:3])   # wet: cool -> warm
tmin,tmax=sm.MAT_C.min(),sm.MAT_C.max(); pmin,pmax=sm.MAP_mm.min(),sm.MAP_mm.max()
def biv_color(mat,map_):
    x=(mat-tmin)/(tmax-tmin); y=(map_-pmin)/(pmax-pmin)
    return tuple((1-x)*(1-y)*_c00+x*(1-y)*_c10+(1-x)*y*_c01+x*y*_c11)
sm["biv"]=[biv_color(r.MAT_C,r.MAP_mm) for r in sm.itertuples()]
p48=p48.merge(sm[["siteID","biv","labcol"]],on="siteID"); pak=pak.merge(sm[["siteID","biv","labcol"]],on="siteID")
# displaced symbol positions for clustered sites (metres), with connectors to true locations
DOTOFF={"ABBY":(-80e3,70e3),"WREF":(80e3,-70e3),"SOAP":(-80e3,-60e3),"TEAK":(80e3,60e3),
        "UNDE":(0,120e3),"STEI":(-120e3,-80e3),"TREE":(120e3,-80e3),
        "BLAN":(-90e3,80e3),"SCBI":(-90e3,-80e3),"SERC":(100e3,0),
        "TALL":(70e3,90e3),"DELA":(-100e3,10e3),"LENO":(40e3,-100e3),
        "ORNL":(-80e3,70e3),"GRSM":(90e3,-60e3)}
AKOFF={"HEAL":(-150e3,0),"BONA":(110e3,100e3),"DEJU":(110e3,-130e3)}
def draw_sites(axh,dd,offmap,ssize):
    xs=[];ys=[]
    for _,r in dd.iterrows():
        dx,dy=offmap.get(r.siteID,(0,0)); x_,y_=r.geometry.x+dx,r.geometry.y+dy
        if dx or dy: axh.plot([r.geometry.x,x_],[r.geometry.y,y_],color="#9a9a9a",lw=0.5,zorder=2)
        xs.append(x_); ys.append(y_)
    axh.scatter(xs,ys,c=list(dd.biv),s=ssize,edgecolor="#222222",linewidth=0.6,zorder=3)
    return xs,ys
fig=plt.figure(figsize=(W*0.82,3.35))
axu=fig.add_axes([0.235,0.01,0.755,0.97])   # CONUS main
axk=fig.add_axes([0.010,0.55,0.205,0.42])   # Alaska, top-left
# --- Alaska panel (own Albers) ---
ak_dom.plot(ax=axk,color="#f4f4f2",edgecolor="none")
ak_dom.boundary.plot(ax=axk,color="#b3b3ad",lw=0.4)
akx,aky=draw_sites(axk,pak,AKOFF,30)
autolab(axk,akx,aky,pak.siteID.values,fs=4.8,colors=pak.labcol.values)
axk.set_xlim(-0.9e6,1.7e6); axk.set_ylim(0.35e6,2.45e6)
axk.set_xticks([]); axk.set_yticks([])
for sp in axk.spines.values(): sp.set_edgecolor("#bbbbbb"); sp.set_linewidth(0.6)
axk.text(0.5,1.03,"Alaska (D18–D19)",transform=axk.transAxes,ha="center",va="bottom",fontsize=6.5,color="#555555")
# --- CONUS panel ---
conus_dom.plot(ax=axu,color="#f4f4f2",edgecolor="none")
conus_dom.boundary.plot(ax=axu,color="#b3b3ad",lw=0.5)
for _,r in conus_dom.iterrows():
    pt=r.geometry.representative_point()
    axu.text(pt.x,pt.y,f"D{int(r.DomainID):02d}",fontsize=4.3,color="#a8a8a2",ha="center",va="center",zorder=1)
usx,usy=draw_sites(axu,p48,DOTOFF,34)
autolab(axu,usx,usy,p48.siteID.values,fs=5.0,colors=p48.labcol.values)
axu.set_axis_off()
# label-colour legend (forest type): left column, between the Alaska panel and the climate gradient
fig.text(0.055,0.525,"Site label = forest type",fontsize=6.2,color="#666666")
for i,(ft,c) in enumerate([("Needleleaf",FTC["Needleleaf"]),("Mixed",FTC["Mixed"]),("Broadleaf",FTC["Broadleaf"])]):
    fig.text(0.062,0.487-0.038*i,ft,fontsize=6.2,fontweight="bold",color=c)
# --- continuous bivariate legend (MAT x MAP), bottom-left of the figure (below Alaska) ---
_figw,_figh=fig.get_size_inches(); _lh=0.30; _lw=_lh*_figh/_figw   # square legend (in inches)
lax=fig.add_axes([0.050,0.09,_lw,_lh])
NG=200
grid=np.zeros((NG,NG,3))
for j,y in enumerate(np.linspace(0,1,NG)):
    for i,x in enumerate(np.linspace(0,1,NG)):
        grid[j,i]=(1-x)*(1-y)*_c00+x*(1-y)*_c10+(1-x)*y*_c01+x*y*_c11
lax.imshow(grid,origin="lower",aspect="auto",extent=[tmin,tmax,pmin,pmax])
lax.set_xticks([0,10,20]); lax.set_yticks([500,1500,2500])
lax.tick_params(labelsize=5,length=1.5,pad=1)
lax.set_xlabel("MAT (°C)",fontsize=6,labelpad=1); lax.set_ylabel("MAP (mm)",fontsize=6,labelpad=1)
for sp in lax.spines.values(): sp.set_visible(False)
fig.savefig(FIG+"/Figure_2.png",bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig2: sites={len(sm)} (26), MAT {sm.MAT_C.min():.1f}~{sm.MAT_C.max():.1f}C, MAP {sm.MAP_mm.min():.0f}~{sm.MAP_mm.max():.0f}mm")

# ================= FIGURE 3 : plot-level framework =================
fig=plt.figure(figsize=(W,9.4)); gs=fig.add_gridspec(3,2,hspace=0.46,wspace=0.34,height_ratios=[1.05,1,1])
x=np.arange(4); w=0.58
# (a) sequential
ax=fig.add_subplot(gs[0,0])
state=[nm.loc[r,"R2_M1_beyond"] for r in RL]; dyn=[nm.loc[r,"R2_M2_beyond"]-nm.loc[r,"R2_M1_beyond"] for r in RL]
ax.bar(x,state,w,label="Canopy state",color=C_STR)
ax.bar(x,dyn,w,bottom=state,label="+ Structural dynamics",color=C_DYN)
for i,r in enumerate(RL):
    if state[i]>0.02: ax.text(i,state[i]/2,f"{state[i]:.3f}",ha="center",va="center",color="white",fontsize=6.6)
    ax.text(i,state[i]+dyn[i]+0.004,f"+{dyn[i]:.3f} {star(wb.loc[r,'p_dyn_wildboot'])}",ha="center",va="bottom",color=C_DYN,fontsize=6.6)
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylabel("R² beyond domain (sequential)")
ax.set_ylim(0,0.185); ax.legend(frameon=False,loc="upper right")
panel(ax,"(a)",-0.16,1.00)
# (b) partition
ax=fig.add_subplot(gs[0,1])
blocks=[("unique_structure","Structure",C_STR),("unique_spectral","Spectral",C_SPE),("unique_dynamics","Dynamics",C_DYN),("shared_RS","Shared",C_SHA)]
for i,r in enumerate(RL):
    b=0
    for col,lab,c in blocks:
        val=vp.loc[r,col]; ax.bar(i,val,w,bottom=b,color=c,label=lab if i==0 else None)
        if val>=0.010:
            txt=f"{val:.3f}"+(f" {star(ps.loc[r,col])}" if col in ps.columns else "")
            ax.text(i,b+val/2,txt,ha="center",va="center",fontsize=6.2,color="white" if c in(C_STR,C_DYN) else "#333333")
        b+=val
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylabel("Unique R² beyond domain (semi-partial)")
ax.set_ylim(0,0.165); ax.legend(frameon=False,loc="upper right")
panel(ax,"(b)",-0.16,1.00)
# (c,d) forests
def forest(sub_ax,resp_label):
    rr=BAYES_RESP[resp_label]; sub=bc[bc.response==rr].set_index("predictor"); order=list(reversed(keep)); ys=np.arange(len(order))
    for yi,pred in zip(ys,order):
        row=sub.loc[pred]; b,lo,hi=row["beta"],row["hdi_2.5"],row["hdi_97.5"]; cred=bool(row["credible"]); c=CBLK[blk(pred)]
        sub_ax.plot([lo,hi],[yi,yi],color=c,lw=2.0 if cred else 1.0,alpha=1.0 if cred else 0.35,solid_capstyle="round",zorder=2)
        sub_ax.scatter([b],[yi],s=26 if cred else 12,color=c,edgecolor="black" if cred else "none",linewidth=0.6,alpha=1.0 if cred else 0.4,zorder=3)
        if cred: sub_ax.text(hi+0.02 if b>0 else lo-0.02,yi,f"{b:+.2f}",va="center",ha="left" if b>0 else "right",fontsize=6.2,color=c)
    sub_ax.axvline(0,color="#999999",lw=0.7,ls="--",zorder=1)
    sub_ax.set_yticks(ys); sub_ax.set_yticklabels([NICE[p] for p in order],fontsize=6.8)
    for tl,p in zip(sub_ax.get_yticklabels(),order): tl.set_color(CBLK[blk(p)])
    sub_ax.set_xlabel("Standardized β (95% HDI)"); sub_ax.set_xlim(-0.55,0.66)
    m=bm.loc[rr]
    sub_ax.set_title(f"{resp_label}   (n = {int(m['n'])})",fontsize=7.6,loc="left",color="#444444")
FOREST_ORDER=[("Hill q1","Hill_q1","(c)"),("Hill q2","Hill_q2","(d)"),
              ("Turnover","LCBD_turnover_rare","(e)"),("Nestedness","LCBD_nestedness_rare","(f)")]
BAYES_RESP={lbl:rr for lbl,rr,_ in FOREST_ORDER if rr in set(bc.response)}
slots=[gs[1,0],gs[1,1],gs[2,0],gs[2,1]]
for (lbl,rr,letter),slot in zip(FOREST_ORDER,slots):
    if rr not in set(bc.response): continue
    axf=fig.add_subplot(slot); forest(axf,lbl); panel(axf,letter,-0.30,1.02)
fig.savefig(FIG+"/Figure_3.png",bbox_inches="tight"); plt.close(fig)
seq=[round(d,3) for d in dyn]
print(f"VERIFY Fig3: seq dyn={seq} (text 0.042/0.039/0.013/0.056) | unique str q1={vp.loc['Hill q1','unique_structure']:.3f} (0.081), dyn nest={vp.loc['LCBD nestedness','unique_dynamics']:.3f} (0.056)")
cred_n=bc[(bc.response=='LCBD_nestedness_rare')&(bc.credible)&(bc.block=='dynamics')]
print(f"VERIFY Fig3: Bayes nest credible dyn = {[(r.predictor,r.beta) for r in cred_n.itertuples()]}")

# ================= FIGURE 4 : variation hypothesis =================
DCACHE=os.path.join(RES,"_fig4_dist_cache.npz")
if os.path.exists(DCACHE):
    _z=np.load(DCACHE); y,xst,xse,nsite=_z["y"],_z["xst"],_z["xse"],int(_z["nsite"])
else:
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
    for s_,g in fz.groupby('siteID'):
        plots=g.plotID.tolist(); sad=dfv[(dfv.siteID==s_)&(dfv.plotID.isin(plots))].groupby(['plotID','taxonID']).size().unstack(fill_value=0)
        common=[p for p in plots if p in sad.index]; gg=g.set_index('plotID').loc[common] if common else None
        if not common or len(common)<6 or gg[['z_'+c for c in Scols]].isna().any().any() or gg[['z_'+c for c in Pcols]].isna().any().any(): continue
        nsite+=1
        Y.append(pdist(sad.loc[common].values,metric='braycurtis'))
        XSt.append(pdist(gg[['z_'+c for c in Scols]].values)); XSe.append(pdist(gg[['z_'+c for c in Pcols]].values))
    y=np.concatenate(Y); xst=np.concatenate(XSt); xse=np.concatenate(XSe)
    np.savez_compressed(DCACHE,y=y,xst=xst,xse=xse,nsite=nsite)
rst=np.corrcoef(rankdata(xst),rankdata(y))[0,1]; rse=np.corrcoef(rankdata(xse),rankdata(y))[0,1]
fig,ax=plt.subplots(1,2,figsize=(W,3.55),sharey=True)
hbs=[]
for a,(xx,r,xlab,letter) in zip(ax,[(xst,rst,"Structural distance (Euclidean)","(a)"),(xse,rse,"Spectral distance (Euclidean)","(b)")]):
    xmax=float(np.quantile(xx,0.995))
    hb=a.hexbin(xx,y,gridsize=(34,26),extent=(0,xmax,0,1.0),cmap="Greys",bins="log",mincnt=1,linewidths=0)
    hbs.append(hb)
    bins=np.quantile(xx,np.linspace(0,1,11)); idx=np.digitize(xx,bins)
    bx=[xx[idx==k].mean() for k in range(1,11) if (idx==k).sum()>=50]; by=[y[idx==k].mean() for k in range(1,11) if (idx==k).sum()>=50]
    a.plot(bx,by,'o-',color=C_ACC,ms=4.2,lw=1.6,zorder=4,label="Decile-bin mean")
    a.set_xlabel(xlab); a.set_xlim(0,xmax); a.set_ylim(-0.02,1.04)
    a.text(0.03,0.045,f"Mantel r = {r:+.2f} (Spearman), p = 0.001",transform=a.transAxes,fontsize=7.2,color="#222222")
    panel(a,letter,-0.14 if letter=="(a)" else -0.06,1.00)
ax[0].set_ylabel("Compositional dissimilarity (Bray–Curtis)")
ax[0].legend(frameon=False,loc="upper right",fontsize=7,bbox_to_anchor=(1.0,0.93))
cb=fig.colorbar(hbs[1],ax=ax,shrink=0.85,pad=0.015,aspect=26)
cb.set_label("Number of plot pairs (log scale)",fontsize=7.5); cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
fig.savefig(FIG+"/Figure_4.png",bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig4: Spearman struct={rst:+.2f} (+0.32), spectral={rse:+.2f} (+0.30), pairs={len(y)} (8,813), sites={nsite} (23)")

# ================= FIGURE 5 : species-energy =================
g3=pd.read_csv(D+"/plot_pml_gpp_ts_26.csv")[["plotID","pml_gpp"]]
gm=pd.read_csv(D+"/plot_modis_gpp_26.csv")[["plotID","modis_gpp"]]
site=(f.merge(g3,on="plotID").merge(gm,on="plotID",how="left")
       .groupby("siteID").agg(Hill_q1=("Hill_q1","mean"),pml=("pml_gpp","mean"),modis=("modis_gpp","mean")).reset_index())
tw=pd.read_csv(D+"/site_tower_gpp_26.csv")[["siteID","tower_gpp"]]; site=site.merge(tw,on="siteID",how="left")
tv=site.dropna(subset=["tower_gpp"])
rt,pt=st.pearsonr(tv.tower_gpp,tv.Hill_q1)
t2=tv.assign(z=(tv.tower_gpp-tv.tower_gpp.mean())/tv.tower_gpp.std()); qt=smf.ols("Hill_q1 ~ z + I(z**2)",t2).fit().pvalues["I(z ** 2)"]
rp_,pp_=st.pearsonr(site.pml,site.Hill_q1); rvt=st.pearsonr(tv.pml,tv.tower_gpp)[0]
sm_=site.dropna(subset=["modis"]); rm_,pm_=st.pearsonr(sm_.modis,sm_.Hill_q1); rmt=st.pearsonr(tv.modis,tv.tower_gpp)[0]
fig,ax=plt.subplots(1,2,figsize=(W*0.72,2.9),sharey=True)
for a,(d3,xc,xl,note,col,letter) in zip(ax,[
    (tv,"tower_gpp","Tower GPP (gC m$^{-2}$ yr$^{-1}$)",f"r = {rt:+.2f}, p = {pt:.3f}\nquadratic ns (p = {qt:.2f})",C_PRO,"(a)"),
    (site,"pml","Site-mean PML-V2 GPP",f"r = {rp_:+.2f}, p = {pp_:.3f}\nvs tower r = +{rvt:.2f}","#2a6f4e","(b)")
    ]):
    a.scatter(d3[xc],d3.Hill_q1,s=24,color=col,edgecolor="white",linewidth=0.5,zorder=3)
    for _,row in d3.iterrows(): a.annotate(row.siteID,(row[xc],row.Hill_q1),fontsize=4.0,color="#666666",xytext=(2.2,1.8),textcoords="offset points")
    b1,b0=np.polyfit(d3[xc],d3.Hill_q1,1); xr=np.linspace(d3[xc].min(),d3[xc].max(),50); a.plot(xr,b0+b1*xr,color="#222222",lw=1.3)
    a.set_xlabel(xl,fontsize=7.8)
    a.text(0.03,0.86,note,transform=a.transAxes,fontsize=6.8,color="#222222")
    panel(a,letter,-0.20 if letter=="(a)" else -0.08,1.00)
ax[0].set_ylabel("Site-mean tree diversity (Hill q1)")
fig.tight_layout(); fig.savefig(FIG+"/Figure_5.png",bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig5: tower r={rt:+.2f} p={pt:.3f} (+0.59,0.004), quad p={qt:.2f} (0.87), PML r={rp_:+.2f} (+0.50), PML-tower +{rvt:.2f} (+0.86), n_tower={len(tv)} (22) | MODIS supplement-only r={rm_:+.2f}, vs tower +{rmt:.2f}")

# ================= FIGURE 6 : context moderation =================
# (a) all FDR-surviving interaction terms from the full screen (168 terms; incl. spectral x land use)
o4=pd.read_csv(RES+"/o4_interactions_pml.csv"); sig=o4[o4.q<0.05].copy()
zcrit=st.norm.ppf(1-sig.p/2); sig["se_appx"]=np.abs(sig.beta_int)/zcrit    # Wald-reconstructed SE
RSN={"LAI_mean":"LAI","VCI_mean":"VCI","EVI_mean":"EVI","pml_gpp":"PML GPP"}
CTXN={"severity":"severity","recency":"recency","lc_shannon":"land-cover diversity","lc_edge":"edge density",
      "forest_frac":"forest fraction","lc_forest_frac":"forest fraction","stand_age_gami":"stand age"}
RESPN={"Hill q1":"Hill q1","Hill q2":"Hill q2","turnover":"Turnover","nestedness":"Nestedness","Turnover":"Turnover","Nestedness":"Nestedness"}
FAMC={"disturbance":C_STR,"land use":"#7d5ba6","stand age":C_PRO}
sig["label"]=[f"{RSN.get(r.rs,r.rs)} × {CTXN.get(r.context,r.context)} → {RESPN.get(r.response,r.response)}" for r in sig.itertuples()]
sig=sig.sort_values(["family","beta_int"]).reset_index(drop=True)
ss=pd.read_csv(RES+"/simple_slopes_pml.csv"); z=[-1,0,1]
gA=ss[(ss.rp=="Hill q1")&(ss.rs=="pml_gpp")&(ss.ctx=="stand_age_gami")].sort_values("ctx_z")
gB=ss[(ss.rp=="Hill q1")&(ss.rs=="VCI_mean")&(ss.ctx=="severity")].sort_values("ctx_z")
fig,ax=plt.subplots(1,3,figsize=(W*1.15,3.0),gridspec_kw={"width_ratios":[1.35,1,1],"wspace":0.32})
a=ax[0]; ysig=np.arange(len(sig))
for yi,r in zip(ysig,sig.itertuples()):
    c=FAMC.get(r.family,"#444444")
    a.plot([r.beta_int-1.96*r.se_appx,r.beta_int+1.96*r.se_appx],[yi,yi],color=c,lw=1.6,solid_capstyle="round")
    a.scatter([r.beta_int],[yi],s=22,color=c,edgecolor="black",linewidth=0.5,zorder=3)
    qtxt="q < 0.001" if r.q<0.001 else f"q = {r.q:.3f}"
    a.text(0.995,yi,qtxt,transform=a.get_yaxis_transform(),fontsize=5.8,color="#777777",va="center",ha="right")
a.axvline(0,color="#999999",lw=0.7,ls="--")
a.set_xlim(-0.245,0.175)
a.set_yticks(ysig); a.set_yticklabels(sig.label,fontsize=6.4)
a.set_xlabel("Interaction β (±95% CI); q = FDR-adjusted p")
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([0],[0],color=FAMC["disturbance"],lw=2,label="Disturbance"),
                  Line2D([0],[0],color=FAMC["land use"],lw=2,label="Land use")],
         frameon=False,loc="upper left",fontsize=6.5,handlelength=1.4,borderaxespad=0.1)
panel(a,"(a)",-0.52,1.02)
# (b,c) simple slopes, separate panels (aspect matched to the previous combined panel)
for a,(g4,col,xt,note,letter) in zip(ax[1:],[
    (gA,C_PRO,["−1 SD\n(younger)","Mean","+1 SD\n(older)"],"MODIS GPP → Hill q1\nby stand age (ns, q = 0.26)","(b)"),
    (gB,C_STR,["−1 SD\n(least dist.)","Mean","+1 SD\n(most dist.)"],"VCI → Hill q1\nby disturbance severity","(c)")]):
    a.errorbar(z,g4.slope,yerr=1.96*g4.se,marker="o",capsize=3,color=col,lw=1.5,ms=5)
    tops=(g4.slope+1.96*g4.se).values; bots=(g4.slope-1.96*g4.se).values
    for zz,sl,tp,pv in zip(z,g4.slope,tops,g4.p):
        a.annotate(f"{sl:+.2f}{'*' if pv<0.05 else ''}",(zz,tp),textcoords="offset points",xytext=(0,3),color=col,fontsize=6.6,ha="center")
    a.axhline(0,color="#999999",lw=0.7,ls="--")
    a.set_xticks(z); a.set_xticklabels(xt,fontsize=6.6); a.set_xlim(-1.55,1.55)
    a.set_ylim(min(-0.02,bots.min()-0.04),tops.max()+0.16)
    a.text(0.04,0.975,note,transform=a.transAxes,fontsize=6.8,color="#444444",va="top")
    panel(a,letter,-0.20,1.02)
ax[1].set_ylabel("Simple slope on Hill q1 (per 1 SD)",fontsize=8)
fig.tight_layout(); fig.savefig(FIG+"/Figure_6.png",bbox_inches="tight"); plt.close(fig)
print(f"VERIFY Fig6: FDR-sig interactions={len(sig)} (5: 2 disturbance + 3 land use) | MODISxAge slopes={[round(v,2) for v in gA.slope]} | VCIxSev={[round(v,2) for v in gB.slope]}")
print("ALL FIGURES SAVED ->",FIG)
