# -*- coding: utf-8 -*-
"""
Recompute the plot-level GPP suppression diagnostic (Results 3.5) on the current
26-site sample with a SINGLE consistent scale (fully standardized: both predictors
and response z-scored), because the previous table mixed raw-scale "alone" models
with standardized "joint" models. Domain FE + site-clustered SE, same retained
predictors as the plot-level framework.
Out: results/O0_framework/gpp_suppression_26.csv
"""
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D, RES = os.path.join(BASE,"data"), os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
S  = ["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]; P = ["EVI_mean"]
Dn = ["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]
G  = ["modis_gpp","pml_gpp"]

f = pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv")); f = f[f.sample_coverage>=0.9].copy()
f = f.merge(pd.read_csv(os.path.join(D,"plot_modis_gpp_26.csv"))[["plotID","modis_gpp"]],on="plotID",how="left")
f = f.merge(pd.read_csv(os.path.join(D,"plot_pml_gpp_ts_26.csv"))[["plotID","pml_gpp"]],on="plotID",how="left")

rows=[]
for resp in ["LCBD_turnover_rare","Hill_q1"]:
    d = f[[resp,"domain","siteID"]+S+P+Dn+G].dropna().copy()
    for c in S+P+Dn+G: d[c]=(d[c]-d[c].mean())/d[c].std()
    y = np.log(d[resp]) if resp.startswith("Hill") else d[resp].astype(float)
    d["y"]=(y-y.mean())/y.std()
    r_gpp = d.modis_gpp.corr(d.pml_gpp)
    for lab,pr in [("modis_only",["modis_gpp"]),("pml_only",["pml_gpp"]),("joint",G)]:
        m = smf.ols("y ~ C(domain) + "+" + ".join(S+P+Dn+pr), d).fit(
                cov_type="cluster", cov_kwds={"groups": d.siteID})
        for g in pr:
            rows.append({"response":resp,"model":lab,"product":g,"n":len(d),
                         "beta_fullstd":round(float(m.params[g]),4),
                         "p_clustered":round(float(m.pvalues[g]),4),
                         "r_modis_pml":round(float(r_gpp),4)})
out = pd.DataFrame(rows)
out.to_csv(os.path.join(RES,"gpp_suppression_26.csv"), index=False)
print(out.to_string(index=False))
print("\nsaved -> gpp_suppression_26.csv")
