import os, warnings, numpy as np, pandas as pd, scipy.stats as st; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
RES=os.path.join(BASE,"papers","UNIFIED","results","O0_framework"); D=os.path.join(BASE,"data"); FIG=os.path.join(BASE,"papers","UNIFIED","figures")
vp=pd.read_csv(RES+"/variance_partition_sd.csv").set_index("response")
nm=pd.read_csv(RES+"/nested_models_sd.csv").set_index("response")
df=pd.read_csv(D+"/FINAL_v2_pooled_26.csv"); df=df[df.sample_coverage>=0.9].copy()
ns=df.siteID.nunique()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];keep=S+P+Dn
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"};LOG={"Hill_q1","Hill_q2"}
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<5e-2 else "n.s."
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
# block-level wildboot p for structure / spectral / dynamics (each tested in full M2)
pstar={}
for rp,lab in RESP.items():
    d=df[[rp,"domain","siteID"]+keep].dropna().copy()
    for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
    y=np.log(d[rp]) if rp in LOG else d[rp].astype(float); d["y"]=(y-y.mean())/y.std()
    pstar[lab]={"unique_structure":wcb(d,["C(domain)"]+P+Dn,S),
                "unique_spectral":wcb(d,["C(domain)"]+S+Dn,P),
                "unique_dynamics":wcb(d,["C(domain)"]+S+P,Dn)}
    print(f"  {lab}: str p={pstar[lab]['unique_structure']:.3f}, spe p={pstar[lab]['unique_spectral']:.3f}, dyn p={pstar[lab]['unique_dynamics']:.3f}",flush=True)

# ---------- FIG 1 ----------
RL=["Hill q1","Hill q2","LCBD turnover","LCBD nestedness"];LAB=["Hill q1\n(alpha)","Hill q2\n(alpha)","LCBD\nturnover","LCBD\nnestedness"]
fig,ax=plt.subplots(figsize=(9.5,6.3)); x=np.arange(4); w=0.62
blocks=[("unique_structure","structure","#08519c"),("unique_spectral","spectral","#9ecae1"),("unique_dynamics","dynamics","#41ab5d"),("shared_RS","shared","#bdbdbd")]
for i,r in enumerate(RL):
    b=0
    for col,lab,c in blocks:
        val=vp.loc[r,col]; ax.bar(i,val,w,bottom=b,color=c,label=lab if i==0 else None)
        if val>0.004:
            txt=f"{val:.3f}"+(f" {star(pstar[r][col])}" if col in pstar[r] else "")
            ax.text(i,b+val/2,txt,ha="center",va="center",fontsize=8,color="white" if col in("unique_structure","unique_dynamics") else "black")
        b+=val
ax.set_xticks(x);ax.set_xticklabels(LAB,fontsize=11);ax.set_ylabel("Unique R² beyond domain (semi-partial)",fontsize=11)
na=int(nm.loc["Hill q1","n"]);nb=int(nm.loc["LCBD turnover","n"])
ax.set_title(f"Plot-level: unique contribution of each remote-sensing dimension\nstructure → alpha; dynamics → richness difference (nestedness); turnover weak\n(n={na} alpha / {nb} beta, {ns} sites; stars = wild cluster bootstrap block tests)",fontweight="bold",fontsize=11)
ax.legend(fontsize=10,title="dimension");fig.tight_layout();fig.savefig(FIG+"/F1_plot_partition.png",dpi=200,bbox_inches="tight");print("saved F1")

# ---------- FIG 2: PML site-mean (26 sites) x-axis, tower validation ----------
g=pd.read_csv(D+"/plot_pml_gpp_ts_26.csv")[["plotID","pml_gpp"]]
site=df.merge(g,on="plotID").groupby("siteID").agg(Hill_q1=("Hill_q1","mean"),pml=("pml_gpp","mean")).reset_index()
tw=pd.read_csv(D+"/site_tower_gpp_26.csv")[["siteID","tower_gpp"]]; site=site.merge(tw,on="siteID",how="left")
r,p=st.pearsonr(site.pml,site.Hill_q1)
s2=site.assign(z=(site.pml-site.pml.mean())/site.pml.std()); quadp=smf.ols("Hill_q1 ~ z + I(z**2)",s2).fit().pvalues["I(z ** 2)"]
rt=st.pearsonr(site.dropna(subset=['tower_gpp']).pml,site.dropna(subset=['tower_gpp']).tower_gpp)[0]
fig,ax=plt.subplots(figsize=(7.6,6))
ax.scatter(site.pml,site.Hill_q1,s=95,color="#00695c",edgecolor="white",linewidth=0.8,zorder=3)
for _,row in site.iterrows(): ax.annotate(row.siteID,(row.pml,row.Hill_q1),fontsize=6.5,alpha=0.6,xytext=(4,3),textcoords="offset points")
b1,b0=np.polyfit(site.pml,site.Hill_q1,1); xr=np.linspace(site.pml.min(),site.pml.max(),50); ax.plot(xr,b0+b1*xr,color="black",lw=2)
ax.set_xlabel("Site-mean PML-V2 GPP (2000–2024)",fontsize=11);ax.set_ylabel("Site-mean tree diversity (Hill q1)",fontsize=11)
ax.set_title(f"Site-level species–energy (PML-V2, all {len(site)} sites)\nr = {r:+.2f}, p = {p:.3f} — monotonic (quadratic n.s., p = {quadp:.2f})\n[PML validated vs eddy-covariance tower: r = +{rt:.2f}, n=22]",fontweight="bold",fontsize=11)
fig.tight_layout();fig.savefig(FIG+"/F2_site_species_energy.png",dpi=200,bbox_inches="tight");print("saved F2")
