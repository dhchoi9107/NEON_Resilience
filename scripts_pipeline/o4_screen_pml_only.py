# -*- coding: utf-8 -*-
"""
O4 context-interaction screen using ONE satellite GPP product (PML-V2).
MODIS MOD17 is dropped from the framework because its light-use-efficiency
formulation is driven by fPAR (a greenness quantity), which conflicts with the
paper's argument for greenness-independent productivity; PML-V2 is a
Penman-Monteith-Leuning model calibrated against eddy-covariance data.
Dropping a predictor changes the multiple-testing burden, so the BH correction
is recomputed on the reduced screen (144 terms instead of 168).
Also recomputes the simple slopes reported in Results 3.6.
Out: o4_interactions_pml.csv, simple_slopes_pml.csv
"""
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D    = os.path.join(BASE, "data")
OUT  = os.path.join(BASE, "papers", "UNIFIED", "results", "O0_framework")

df = pd.read_csv(os.path.join(D, "FINAL_v2_pooled_26.csv"))
df = df.merge(pd.read_csv(os.path.join(D, "plot_pml_gpp_ts_26.csv"))[["plotID","pml_gpp"]], on="plotID", how="left")
df = df.merge(pd.read_csv(os.path.join(D, "plot_stand_age_gami_26.csv"))[["plotID","stand_age_gami"]], on="plotID", how="left")
df = df.merge(pd.read_csv(os.path.join(D, "plot_disturbance_robust_26.csv"))[["plotID","severity","recency"]], on="plotID", how="left")
df = df.merge(pd.read_csv(os.path.join(D, "plot_landuse_het_26.csv"))[["plotID","lc_edge","lc_shannon","lc_forest_frac"]], on="plotID", how="left")
df = df[df.sample_coverage >= 0.90].copy()

RS   = ["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean","EVI_mean","pml_gpp"]
RESP = {"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"turnover","LCBD_nestedness_rare":"nestedness"}
LOG  = {"Hill_q1","Hill_q2"}
FAM  = {"stand age":["stand_age_gami"], "disturbance":["severity","recency"],
        "land use":["lc_edge","lc_shannon","lc_forest_frac"]}
def z(s): return (s - s.mean())/s.std()

rows = []
for fam, ctxs in FAM.items():
    for ctx in ctxs:
        for rp, lab in RESP.items():
            for rs in RS:
                d = df[[rp,"domain","siteID",rs,ctx]].dropna().copy()
                if d.siteID.nunique() < 5 or len(d) < 40: continue
                d["y"]  = np.log(d[rp]) if rp in LOG else d[rp].astype(float)
                d["zR"], d["zC"] = z(d[rs]), z(d[ctx])
                try:
                    m = smf.ols("y ~ C(domain) + zR + zC + zR:zC", d).fit(
                            cov_type="cluster", cov_kwds={"groups": d["siteID"]})
                    rows.append(dict(family=fam, context=ctx, response=lab, rs=rs,
                                     beta_int=m.params.get("zR:zC", np.nan),
                                     p=m.pvalues.get("zR:zC", np.nan), n=len(d)))
                except Exception: pass
I = pd.DataFrame(rows).dropna(subset=["p"])
I["q"] = np.nan
for fam in FAM:
    msk = I.family == fam
    I.loc[msk, "q"] = multipletests(I.loc[msk, "p"].values, method="fdr_bh")[1]
I["sig"] = I["q"] < 0.05
I.to_csv(os.path.join(OUT, "o4_interactions_pml.csv"), index=False)
print(f"=== O4 screen, PML only: {len(I)} terms, {int(I.sig.sum())} survive BH ===")
print(I[I.sig][["family","context","response","rs","beta_int","p","q"]].round(4).to_string(index=False))

# ---------- simple slopes at -1 SD / mean / +1 SD ----------
def slopes(rp, rs, ctx):
    d = df[[rp,"domain","siteID",rs,ctx]].dropna().copy()
    d["y"] = np.log(d[rp]) if rp in LOG else d[rp].astype(float)
    d["y"] = z(d["y"])          # response standardized too (matches 159_fewcluster_robust_26)
    d["zR"], d["zC"] = z(d[rs]), z(d[ctx])
    m = smf.ols("y ~ C(domain) + zR + zC + zR:zC", d).fit(cov_type="cluster", cov_kwds={"groups": d["siteID"]})
    b, bi = m.params["zR"], m.params.get("zR:zC", 0.0)
    V = m.cov_params()
    out = []
    for zc in (-1, 0, 1):
        est = b + bi*zc
        se  = np.sqrt(V.loc["zR","zR"] + zc**2*V.loc["zR:zC","zR:zC"] + 2*zc*V.loc["zR","zR:zC"])
        from scipy import stats
        out.append(dict(rp=RESP[rp], rs=rs, ctx=ctx, ctx_z=zc, slope=est, se=se,
                        p=2*stats.norm.sf(abs(est/se))))
    return out
sl = []
for rs, ctx in [("pml_gpp","stand_age_gami"), ("VCI_mean","severity"), ("LAI_mean","severity")]:
    sl += slopes("Hill_q1", rs, ctx)
S = pd.DataFrame(sl); S.to_csv(os.path.join(OUT, "simple_slopes_pml.csv"), index=False)
print("\n=== simple slopes (Hill q1) ===")
print(S.round(4).to_string(index=False))
