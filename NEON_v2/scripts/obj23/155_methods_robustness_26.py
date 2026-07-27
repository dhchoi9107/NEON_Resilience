"""
Robustness checks to make revised Methods claims TRUE:
 (1) height collinearity (r) + does adding canopy height improve M3 (ΔAIC, ΔR2)?
 (2) singular mixed model vs clustered-OLS: fixed-effect coefficient agreement
 (3) LCBD bounded-response robustness: beta-regression nested LRT vs OLS conclusions
 (4) count of interaction terms actually tested
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
import patsy
from statsmodels.othermod.betareg import BetaModel

D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"
df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv"))
for f in ["plot_modis_gpp_26.csv","plot_pml_gpp_26.csv"]: df=df.merge(pd.read_csv(os.path.join(D,f)),on="plotID",how="left")
df=df[df.sample_coverage>=0.90].copy()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];G=["modis_gpp","pml_gpp"]
keep=S+P+Dn+G

# ---------- (1) height collinearity + incremental value ----------
print("=== (1) Canopy height collinearity & incremental value ===")
cc=df[["Canopy_Ht_mean","Rugosity_mean","LAI_mean","VCI_mean"]].dropna()
for v in ["Rugosity_mean","LAI_mean","VCI_mean"]:
    print(f"  r(Canopy_Ht, {v}) = {cc['Canopy_Ht_mean'].corr(cc[v]):+.2f}")
d=df[["Hill_q1","domain","siteID"]+keep+["Canopy_Ht_mean"]].dropna().copy()
for c in keep+["Canopy_Ht_mean"]: d[c]=(d[c]-d[c].mean())/d[c].std()
d["y"]=np.log(d["Hill_q1"])
m3 =smf.ols("y ~ C(domain)"+"".join(" + "+c for c in keep),d).fit()
m3h=smf.ols("y ~ C(domain)"+"".join(" + "+c for c in keep+['Canopy_Ht_mean']),d).fit()
print(f"  add Canopy_Ht to M3: dAIC={m3h.aic-m3.aic:+.2f}  dR2={m3h.rsquared-m3.rsquared:+.4f}  p(height)={m3h.pvalues['Canopy_Ht_mean']:.2g}")

# ---------- (2) singular mixed model vs clustered OLS: coefficient agreement ----------
print("\n=== (2) singular mixed (domain fixed + site random) vs OLS+domain: fixed-effect agreement ===")
d2=df[["Hill_q1","domain","siteID"]+keep].dropna().copy()
for c in keep: d2[c]=(d2[c]-d2[c].mean())/d2[c].std()
d2["y"]=np.log(d2["Hill_q1"])
f="y ~ C(domain)"+"".join(" + "+c for c in keep)
ols=smf.ols(f,d2).fit()
mix=smf.mixedlm(f,d2,groups=d2["siteID"]).fit(reml=False,method="lbfgs")
comp=pd.DataFrame({"ols":ols.params,"mixed":mix.fe_params}).dropna()
comp=comp.loc[[i for i in comp.index if i in keep]]  # RS predictors only
print(f"  max |Δcoef| = {(comp.ols-comp.mixed).abs().max():.4f}   corr = {comp.ols.corr(comp.mixed):.4f}   (n={len(comp)} RS predictors)")

# ---------- (3) LCBD: beta-regression nested LRT vs OLS ----------
print("\n=== (3) LCBD bounded-response: beta regression nested LRT (vs OLS) ===")
def squeeze(y):  # Smithson-Verkuilen to open (0,1)
    n=len(y); return (y*(n-1)+0.5)/n
def beta_llf(formula,data):
    y,X=patsy.dmatrices(formula,data,return_type="dataframe")
    return BetaModel(y.iloc[:,0].values, X).fit(disp=0)
def lrt(llf_f,llf_r,k): c=2*(llf_f-llf_r); return c, st.chi2.sf(c,k)
for resp,lab in [("LCBD_turnover_rare","turnover"),("LCBD_nestedness_rare","nestedness")]:
    d3=df[[resp,"domain","siteID"]+keep].dropna().copy()
    for c in keep: d3[c]=(d3[c]-d3[c].mean())/d3[c].std()
    d3["yb"]=squeeze(d3[resp].values)
    f1="yb ~ C(domain)"+"".join(" + "+c for c in S+P)
    f2="yb ~ C(domain)"+"".join(" + "+c for c in S+P+Dn)
    f3="yb ~ C(domain)"+"".join(" + "+c for c in S+P+Dn+G)
    b1,b2,b3=beta_llf(f1,d3),beta_llf(f2,d3),beta_llf(f3,d3)
    c12,p12=lrt(b2.llf,b1.llf,len(Dn)); c23,p23=lrt(b3.llf,b2.llf,len(G))
    print(f"  {lab:10s} beta-reg  M1->M2(dyn) p={p12:.2g}   M2->M3(fun) p={p23:.2g}")

# ---------- (4) interaction term count ----------
print("\n=== (4) interaction terms actually tested (per family) ===")
tot=0
for fn,fam,pcol,tcol in [("obj2_severity_recency.csv","disturbance",None,"term"),
                         ("stand_age_moderation_gami.csv","stand age",None,"rs"),
                         ("obj2_heterogeneity.csv","land use",None,"term")]:
    t=pd.read_csv(os.path.join(R,fn))
    if fn=="obj2_severity_recency.csv": t=t[t.test.str.contains("x_RS|inter|recency",case=False,na=False)]
    if fn=="obj2_heterogeneity.csv": t=t[t.test.str.contains("x|mod|inter",case=False,na=False)] if "test" in t else t.iloc[0:0]
    t=t[~t[tcol].astype(str).str.contains("dhi",case=False,na=False)]
    print(f"  {fam:12s}: {len(t)} interaction terms"); tot+=len(t)
print(f"  TOTAL interaction terms (DHI excluded) = {tot}")
import sys; print("\nPython", sys.version.split()[0], "| statsmodels", __import__("statsmodels").__version__)
