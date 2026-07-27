"""
Climate control (NEON site MAT/MAP). Runs the nested framework under a CLIMATE baseline and
compares with the DOMAIN baseline. Serves both:
  Approach 1 = climate-baseline SENSITIVITY (does RS add beyond climate? esp. GPP/species-energy)
  Approach 2 = climate as the MAIN baseline layer (Climate -> State -> Dynamics -> Function)
OLS + site-clustered SE. Hill log-transformed. Out: papers/UNIFIED/results/O0_framework/climate_*.csv
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor

BASE=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"; D=os.path.join(BASE,"data")
OUT=os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv"))
for f in ["plot_modis_gpp_26.csv","plot_pml_gpp_26.csv"]: df=df.merge(pd.read_csv(os.path.join(D,f)),on="plotID",how="left")
df=df.merge(pd.read_csv(os.path.join(D,"site_climate_neon.csv")),on="siteID",how="left")

S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];G=["modis_gpp","pml_gpp"]
CLIM=["MAT_C","MAP_mm"]                 # energy + water (elevation collinear w/ MAT -> excluded; see VIF)
keep=S+P+Dn+G; df=df[df.sample_coverage>=0.90].copy()
RESP={"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}
LOG={"Hill_q1","Hill_q2"}

# climate VIF (incl elevation for the record)
vc=df[["MAT_C","MAP_mm","elevation_m"]].dropna().copy()
for c in vc: vc[c]=(vc[c]-vc[c].mean())/vc[c].std()
print("climate VIF:", {c:round(variance_inflation_factor(vc.values,i),2) for i,c in enumerate(vc.columns)})
print("corr(MAT,elev)=%.2f\n"%vc["MAT_C"].corr(vc["elevation_m"]))

def lrt(a,b,k): c=2*(a.llf-b.llf); return c,(st.chi2.sf(c,k) if k>0 and c>0 else np.nan)
rows=[]; o3=[]
for rp,lab in RESP.items():
    d=df[[rp,"domain","siteID"]+keep+CLIM].dropna().copy()
    for c in keep+CLIM: d[c]=(d[c]-d[c].mean())/d[c].std()
    d["y"]=np.log(d[rp]) if rp in LOG else d[rp]
    def fit(base,preds): return smf.ols(f"y ~ {base}"+"".join(" + "+c for c in preds),d).fit()
    for bl_name, base in [("domain","C(domain)"), ("climate"," + ".join(CLIM))]:
        m0=fit(base,[]); m1=fit(base,S+P); m2=fit(base,S+P+Dn); m3=fit(base,S+P+Dn+G)
        c12=lrt(m2,m1,len(Dn)); c23=lrt(m3,m2,len(G))
        rows.append(dict(response=lab,baseline=bl_name,n=len(d),
            R2_base=m0.rsquared,R2_M1=m1.rsquared,R2_M2=m2.rsquared,R2_M3=m3.rsquared,
            dR2_dyn=m2.rsquared-m1.rsquared, p_dyn=c12[1], dAIC_dyn=m2.aic-m1.aic,
            dR2_fun=m3.rsquared-m2.rsquared, p_fun=c23[1], dAIC_fun=m3.aic-m2.aic))
    # O3 focus: does GPP add beyond climate ALONE, and is the sign positive? (species-energy beyond climate)
    mc =fit(" + ".join(CLIM),[]); mcg=fit(" + ".join(CLIM),G)
    cg=lrt(mcg,mc,len(G))
    mcg_cl=smf.ols("y ~ "+" + ".join(CLIM)+" + "+" + ".join(G),d).fit(cov_type="cluster",cov_kwds={"groups":d["siteID"]})
    o3.append(dict(response=lab, R2_climate=mc.rsquared, dR2_gpp_beyond_climate=mcg.rsquared-mc.rsquared,
        p_gpp=cg[1], modis_beta=mcg_cl.params["modis_gpp"], modis_p_clustered=mcg_cl.pvalues["modis_gpp"],
        pml_beta=mcg_cl.params["pml_gpp"]))

comp=pd.DataFrame(rows); comp.to_csv(os.path.join(OUT,"climate_nested_compare.csv"),index=False)
o3d=pd.DataFrame(o3); o3d.to_csv(os.path.join(OUT,"climate_species_energy.csv"),index=False)
pd.set_option("display.width",220)
print("=== NESTED under DOMAIN vs CLIMATE baseline (dynamics & function increments) ===")
print(comp[["response","baseline","R2_base","dR2_dyn","p_dyn","dR2_fun","p_fun"]].round(4).to_string(index=False))
print("\n=== O3: does GPP add BEYOND CLIMATE ALONE? (species-energy) ===")
print(o3d.round(4).to_string(index=False))
