# Restructured overall figure: plot-level framework (state + dynamics, NO productivity)
# + site-level species-energy (tower GPP) as productivity at its correct scale.
import os, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"; RES=os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
D=os.path.join(BASE,"data"); OUT=os.path.join(BASE,"papers","UNIFIED","figures","OVERALL_sd.png")
nm=pd.read_csv(RES+"/nested_models_sd.csv").set_index("response")
vp=pd.read_csv(RES+"/variance_partition_sd.csv").set_index("response")
wb=pd.read_csv(RES+"/wildboot_sd.csv").set_index("response")
ns=pd.read_csv(D+"/FINAL_v2_pooled_26.csv").siteID.nunique()
RESP=["Hill q1","Hill q2","LCBD turnover","LCBD nestedness"]; LAB=["Hill q1\n(alpha)","Hill q2\n(alpha)","LCBD\nturnover","LCBD\nnestedness"]
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "n.s."
na=int(nm.loc["Hill q1","n"]); nb=int(nm.loc["LCBD turnover","n"])
fig,ax=plt.subplots(1,3,figsize=(18,5.5))
fig.suptitle(f"Remote sensing & tree diversity — restructured (productivity at site scale)  "
             f"[plot framework n={na} alpha/{nb} beta, {ns} sites; stars=wild cluster bootstrap]",fontweight="bold",fontsize=12.5)
x=np.arange(4); w=0.6
# (A) sequential state -> +dynamics
st_=[nm.loc[r,"R2_M1_beyond"] for r in RESP]; dy=[nm.loc[r,"R2_M2_beyond"]-nm.loc[r,"R2_M1_beyond"] for r in RESP]
ax[0].bar(x,st_,w,label="State (structure+spectral)",color="#2c7fb8")
ax[0].bar(x,dy,w,bottom=st_,label="+ Dynamics (LiDAR trends)",color="#41ab5d")
for i,r in enumerate(RESP):
    ax[0].text(i,st_[i]/2,f"{st_[i]:.3f}",ha="center",va="center",color="white",fontweight="bold",fontsize=9)
    ax[0].text(i,st_[i]+dy[i]/2,f"+{dy[i]:.3f}\n{star(wb.loc[r,'p_dyn_wildboot'])}",ha="center",va="center",color="white",fontsize=8)
ax[0].set_xticks(x);ax[0].set_xticklabels(LAB);ax[0].set_ylabel("Incremental R² beyond domain")
ax[0].set_title("(A) State → +Dynamics (plot level)\nno productivity block",fontweight="bold",fontsize=11);ax[0].legend(fontsize=8.5)
# (B) variance partition
blocks=[("unique_structure","structure","#08519c"),("unique_spectral","spectral","#9ecae1"),("unique_dynamics","dynamics","#41ab5d"),("shared_RS","shared","#bdbdbd")]
bot=np.zeros(4)
for col,lab,c in blocks:
    v=[vp.loc[r,col] for r in RESP]; ax[1].bar(x,v,w,bottom=bot,label=lab,color=c); bot+=np.array(v)
ax[1].set_xticks(x);ax[1].set_xticklabels(LAB);ax[1].set_ylabel("R² beyond domain (partitioned)")
ax[1].set_title("(B) Unique contribution\n(structure→alpha; dynamics→nestedness)",fontweight="bold",fontsize=11);ax[1].legend(fontsize=8.5)
# (C) site-level species-energy: tower GPP vs Hill q1
f=pd.read_csv(D+"/FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9]
sd=f.groupby("siteID").Hill_q1.mean().reset_index().merge(pd.read_csv(D+"/site_tower_gpp_26.csv")[["siteID","tower_gpp"]],on="siteID").dropna()
r,p=st.pearsonr(sd.tower_gpp,sd.Hill_q1)
ax[2].scatter(sd.tower_gpp,sd.Hill_q1,s=70,color="#c62828",edgecolor="white",zorder=3)
b1,b0=np.polyfit(sd.tower_gpp,sd.Hill_q1,1); xr=np.linspace(sd.tower_gpp.min(),sd.tower_gpp.max(),50)
ax[2].plot(xr,b0+b1*xr,color="black",lw=2)
ax[2].set_xlabel("Tower GPP (gC m⁻² yr⁻¹, eddy-covariance)");ax[2].set_ylabel("Site-mean Hill q1")
ax[2].set_title(f"(C) Productivity = site-level species–energy\ntower GPP: r={r:+.2f}, p={p:.3f}, n={len(sd)} (monotonic)",fontweight="bold",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig(OUT,dpi=200,bbox_inches="tight"); print("saved",OUT)
