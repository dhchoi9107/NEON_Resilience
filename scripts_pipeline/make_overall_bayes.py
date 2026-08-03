"""
Overall results figure, Bayesian-reflected (restructured state->dynamics framework, 26 sites).
 (A) Plot-level variance partition: unique R^2 beyond domain per block, wild cluster bootstrap
     block stars. The dynamics->nestedness bar is flagged Bayesian-credible (see B).
 (B) Bayesian multilevel forest — LCBD nestedness: 3 dynamics trends credible (95% HDI excl 0)
     even though the joint dynamics block test is n.s. (the headline Bayesian result).
 (C) Bayesian multilevel forest — Hill q1 alpha: structure dominates.
 (D) Site-level species-energy (PML-V2, tower-validated).
Reads results/O0_framework/{variance_partition_sd,wildboot_sd,bayes_multilevel_coeffs}.csv.
Structure/spectral block wildboot p computed inline (dynamics read from wildboot_sd.csv).
Out: papers/UNIFIED/figures/OVERALL_bayes.png
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
RES=os.path.join(BASE,"papers","UNIFIED","results","O0_framework"); D=os.path.join(BASE,"data"); FIG=os.path.join(BASE,"papers","UNIFIED","figures")
vp=pd.read_csv(RES+"/variance_partition_sd.csv").set_index("response")
wbsd=pd.read_csv(RES+"/wildboot_sd.csv").set_index("response")   # dynamics block p
bc=pd.read_csv(RES+"/bayes_multilevel_coeffs.csv")
df=pd.read_csv(D+"/FINAL_v2_pooled_26.csv"); df=df[df.sample_coverage>=0.9].copy(); ns=df.siteID.nunique()

S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]; P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]; keep=S+P+Dn
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}; LOG={"Hill_q1","Hill_q2"}
RL=["Hill q1","Hill q2","LCBD turnover","LCBD nestedness"]; LAB=["Hill q1\n(alpha)","Hill q2\n(alpha)","LCBD\nturnover","LCBD\nnestedness"]
NICE={"Rugosity_mean":"Rugosity","Vert_CV_mean":"Vert-CV","VCI_mean":"VCI","LAI_mean":"LAI","EVI_mean":"EVI",
      "Rumple_trend":"Rumple trend","Vert_SD_trend":"Vert-SD trend","Vert_CV_trend":"Vert-CV trend","VCI_trend":"VCI trend",
      "FHD_trend":"FHD trend","LAI_trend":"LAI trend","Ht_Ratio_trend":"Ht-ratio trend"}
CBLK={"structure":"#08519c","spectral":"#9ecae1","dynamics":"#41ab5d"}
def blk(c): return "structure" if c in S else "spectral" if c in P else "dynamics"
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "n.s."

# --- inline wild cluster bootstrap block p for structure & spectral (dynamics from CSV) ---
def wcb(d,base_terms,block,B=999,seed=11):
    rng=np.random.RandomState(seed); ff="y ~ "+" + ".join(base_terms+block); fr="y ~ "+" + ".join(base_terms)
    m=smf.ols(ff,d).fit(cov_type="cluster",cov_kwds={"groups":d.siteID})
    Fo=float(m.wald_test(",".join(f"{t}=0" for t in block),use_f=True,scalar=True).statistic)
    r0=smf.ols(fr,d).fit(); fit0,res0=r0.fittedvalues,r0.resid; uniq=np.unique(d.siteID.values); cnt=0
    for _ in range(B):
        w=pd.Series(np.where(rng.rand(len(uniq))<0.5,1,-1),index=uniq)
        dd=d.assign(y=fit0+res0*d.siteID.map(w).values)
        Fb=float(smf.ols(ff,dd).fit(cov_type="cluster",cov_kwds={"groups":dd.siteID}).wald_test(",".join(f"{t}=0" for t in block),use_f=True,scalar=True).statistic)
        if Fb>=Fo: cnt+=1
    return (cnt+1)/(B+1)
CACHE=os.path.join(RES,"_pstar_blocks_cache.csv")
if os.path.exists(CACHE):
    pc=pd.read_csv(CACHE).set_index("response")
    pstar={lab:{"unique_structure":float(pc.loc[lab,"unique_structure"]),"unique_spectral":float(pc.loc[lab,"unique_spectral"]),
                "unique_dynamics":float(pc.loc[lab,"unique_dynamics"])} for lab in RESP.values()}
    print("loaded pstar cache",flush=True)
else:
    pstar={}
    for rp,lab in RESP.items():
        d=df[[rp,"domain","siteID"]+keep].dropna().copy()
        for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
        y=np.log(d[rp]) if rp in LOG else d[rp].astype(float); d["y"]=(y-y.mean())/y.std()
        pstar[lab]={"unique_structure":wcb(d,["C(domain)"]+P+Dn,S),"unique_spectral":wcb(d,["C(domain)"]+S+Dn,P),
                    "unique_dynamics":float(wbsd.loc[lab,"p_dyn_wildboot"])}
        print(f"  {lab}: str={pstar[lab]['unique_structure']:.3f} spe={pstar[lab]['unique_spectral']:.3f} dyn={pstar[lab]['unique_dynamics']:.3f}",flush=True)
    pd.DataFrame([{"response":k,**v} for k,v in pstar.items()]).to_csv(CACHE,index=False)

# ---- which blocks have >=1 Bayesian-credible coefficient (only Hill q1 & nestedness modelled) ----
BAYES_RESP={"Hill q1":"Hill_q1","LCBD nestedness":"LCBD_nestedness_rare"}
bcred={}   # (response_label -> set of blocks with a credible coef)
for lab,rr in BAYES_RESP.items():
    sub=bc[(bc.response==rr)&(bc.credible)]
    bcred[lab]=set(sub.block.unique())

# =========================== FIGURE ===========================
fig=plt.figure(figsize=(15.5,11))
gs=GridSpec(2,2,figure=fig,hspace=0.32,wspace=0.24)
na=int(df.dropna(subset=["Hill_q1"]).shape[0])
fig.suptitle("Remote-sensing dimensions of tree diversity (26 NEON sites) — plot-level framework, Bayesian-reflected",
             fontsize=15,fontweight="bold",y=0.985)

# ---------- (A) variance partition ----------
ax=fig.add_subplot(gs[0,0]); x=np.arange(4); w=0.62
blocks=[("unique_structure","structure","#08519c"),("unique_spectral","spectral","#9ecae1"),("unique_dynamics","dynamics","#41ab5d"),("shared_RS","shared","#bdbdbd")]
for i,r in enumerate(RL):
    b=0
    for col,lab,c in blocks:
        val=vp.loc[r,col]; ax.bar(i,val,w,bottom=b,color=c,label=lab if i==0 else None)
        if val>0.004:
            txt=f"{val:.3f}"+(f" {star(pstar[r][col])}" if col in pstar[r] else "")
            ax.text(i,b+val/2,txt,ha="center",va="center",fontsize=8,color="white" if col in("unique_structure","unique_dynamics") else "black")
        b+=val
# flag the dynamics->nestedness bar as Bayesian-credible
ni=RL.index("LCBD nestedness"); ytop=vp.loc["LCBD nestedness",["unique_structure","unique_spectral","unique_dynamics","shared_RS"]].sum()
ax.annotate("dynamics→nestedness:\nBayesian-credible (panel B),\nblock test n.s.",
            xy=(ni-0.28,vp.loc["LCBD nestedness","unique_structure"]+vp.loc["LCBD nestedness","unique_dynamics"]*0.55),
            xytext=(1.62,0.158),fontsize=8.4,fontweight="bold",color="#1a7a3a",ha="center",va="center",
            arrowprops=dict(arrowstyle="->",color="#1a7a3a",lw=1.5,connectionstyle="arc3,rad=-0.15"))
ax.set_xticks(x); ax.set_xticklabels(LAB,fontsize=10.5); ax.set_ylabel("Unique R² beyond domain (semi-partial)",fontsize=11)
ax.set_ylim(0,ytop+0.11)
ax.set_title("(A) Unique contribution of each dimension\nstars = wild cluster bootstrap joint block test (conservative at few clusters)",fontweight="bold",fontsize=11)
ax.legend(fontsize=9,title="dimension",loc="upper right")

# ---------- forest plots ----------
def forest(ax,resp_label,title):
    rr=BAYES_RESP[resp_label]; sub=bc[bc.response==rr].set_index("predictor")
    order=list(reversed(keep))                     # structure top -> dynamics bottom
    ys=np.arange(len(order))
    for yi,pred in zip(ys,order):
        row=sub.loc[pred]; b,lo,hi=row["beta"],row["hdi_2.5"],row["hdi_97.5"]; cred=bool(row["credible"]); c=CBLK[blk(pred)]
        ax.plot([lo,hi],[yi,yi],color=c,lw=2.6 if cred else 1.4,alpha=1.0 if cred else 0.4,solid_capstyle="round",zorder=2)
        ax.scatter([b],[yi],s=70 if cred else 34,color=c,edgecolor="black" if cred else "none",linewidth=0.8,alpha=1.0 if cred else 0.45,zorder=3)
        if cred: ax.text(hi+0.012 if b>0 else lo-0.012,yi,f"{b:+.2f}",va="center",ha="left" if b>0 else "right",fontsize=7.6,fontweight="bold",color=c)
    ax.axvline(0,color="grey",lw=1,ls="--",zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels([NICE[p] for p in order],fontsize=9)
    ax.set_xlabel("Standardized β (95% HDI)",fontsize=10.5); ax.set_xlim(-0.55,0.65)
    ax.set_title(title,fontweight="bold",fontsize=11)
    # block color y-tick labels
    for tl,p in zip(ax.get_yticklabels(),order): tl.set_color(CBLK[blk(p)])
    m=pd.read_csv(RES+"/bayes_multilevel_meta.csv").set_index("response").loc[rr]
    ax.text(0.5,-0.17,f"0 divergences · R̂ ≤ {m['max_rhat']:.3f} · n = {int(m['n'])} · filled = 95% HDI excludes 0",
            transform=ax.transAxes,ha="center",va="top",fontsize=8,style="italic",color="grey")

axB=fig.add_subplot(gs[0,1]); forest(axB,"LCBD nestedness","(B) Bayesian multilevel — LCBD nestedness\ndynamics trends credible (Ht-ratio, Vert-CV, LAI)")
axD=fig.add_subplot(gs[1,1]); forest(axD,"Hill q1","(D) Bayesian multilevel — Hill q1 (alpha)\nstructure dominates (LAI, VCI) + spectral (EVI)")

# ---------- (C) site-level species-energy ----------
axC=fig.add_subplot(gs[1,0])
g=pd.read_csv(D+"/plot_pml_gpp_ts_26.csv")[["plotID","pml_gpp"]]
site=df.merge(g,on="plotID").groupby("siteID").agg(Hill_q1=("Hill_q1","mean"),pml=("pml_gpp","mean")).reset_index()
tw=pd.read_csv(D+"/site_tower_gpp_26.csv")[["siteID","tower_gpp"]]; site=site.merge(tw,on="siteID",how="left")
r,p=st.pearsonr(site.pml,site.Hill_q1)
s2=site.assign(z=(site.pml-site.pml.mean())/site.pml.std()); quadp=smf.ols("Hill_q1 ~ z + I(z**2)",s2).fit().pvalues["I(z ** 2)"]
sv=site.dropna(subset=['tower_gpp']); rt=st.pearsonr(sv.pml,sv.tower_gpp)[0]
axC.scatter(site.pml,site.Hill_q1,s=85,color="#00695c",edgecolor="white",linewidth=0.8,zorder=3)
for _,row in site.iterrows(): axC.annotate(row.siteID,(row.pml,row.Hill_q1),fontsize=6,alpha=0.55,xytext=(4,3),textcoords="offset points")
b1,b0=np.polyfit(site.pml,site.Hill_q1,1); xr=np.linspace(site.pml.min(),site.pml.max(),50); axC.plot(xr,b0+b1*xr,color="black",lw=2)
axC.set_xlabel("Site-mean PML-V2 GPP (2000–2024)",fontsize=10.5); axC.set_ylabel("Site-mean tree diversity (Hill q1)",fontsize=10.5)
axC.set_title(f"(C) Site-level species–energy (all {len(site)} sites)\nr = {r:+.2f}, p = {p:.3f}; monotonic (quadratic n.s., p={quadp:.2f}); PML vs tower r=+{rt:.2f}",fontweight="bold",fontsize=11)

fig.savefig(FIG+"/OVERALL_bayes.png",dpi=200,bbox_inches="tight"); print("saved",FIG+"/OVERALL_bayes.png")
