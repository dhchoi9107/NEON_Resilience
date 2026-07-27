"""
Unified analytical framework (makes MANUSCRIPT 2.5 true):
 1. VIF screening (< 5) of remote-sensing predictors
 2. Hierarchical nested models  M0 domain -> M1 State -> M2 +Dynamics -> M3 +Function(GPP)
    OLS with NEON domain fixed effects (biogeographic control) + SITE-CLUSTERED robust SEs.
    (A site random intercept was evaluated but its variance collapsed to ~0 once canopy
    predictors were included -> singular; OLS+domain+cluster-SE is the appropriate, finite model.)
    Reported: R2, R2 beyond domain, AIC, dAIC, and likelihood-ratio tests between nested models.
 3. Variance partitioning into unique/shared fractions (structure/spectral/dynamics/productivity)
 4. FDR (Benjamini-Hochberg) on O4 interaction terms (DHI terms excluded)
 5. Moran's I on residuals (kNN weights, permutation) + spatial-covariate robustness refit
Out: papers/UNIFIED/results/O0_framework/*.csv  (+ console summary)
"""
import os, warnings, numpy as np, pandas as pd, scipy.stats as st
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.multitest import multipletests
from sklearn.neighbors import NearestNeighbors

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D = os.path.join(BASE, "data"); R = os.path.join(BASE, "results")
OUT = os.path.join(BASE, "papers", "UNIFIED", "results", "O0_framework"); os.makedirs(OUT, exist_ok=True)
np.random.seed(42)

# ---------- data ----------
df = pd.read_csv(os.path.join(D, "FINAL_v2_pooled_26.csv"))
df = df.merge(pd.read_csv(os.path.join(D, "plot_modis_gpp_26.csv")), on="plotID", how="left")
df = df.merge(pd.read_csv(os.path.join(D, "plot_pml_gpp_26.csv")), on="plotID", how="left")
df = df.merge(pd.read_csv(os.path.join(D, "plot_lonlat_26.csv"))[["plotID", "lon", "lat"]], on="plotID", how="left")
df = df[df.sample_coverage >= 0.90].copy()

STRUCT = [f"{m}_mean" for m in ["Canopy_Ht","Max_Ht","Rumple","Rugosity","Deep_Gap","Vert_SD",
          "Vert_CV","Gini","VCI","FHD","LAI","Q95","Ht_Ratio"]]
SPEC   = [f"{m}_mean" for m in ["NDVI","EVI","ARVI","SAVI"]]
DYN    = [f"{m}_trend" for m in ["Canopy_Ht","Max_Ht","Rumple","Rugosity","Deep_Gap","Vert_SD",
          "Vert_CV","Gini","VCI","FHD","LAI","Q95","Ht_Ratio"]]   # STRUCTURAL trends only
PROD   = ["modis_gpp", "pml_gpp"]
RESP   = {"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}
ALLP = STRUCT + SPEC + DYN + PROD
blk = {**{c:"structure" for c in STRUCT}, **{c:"spectral" for c in SPEC},
       **{c:"dynamics" for c in DYN}, **{c:"productivity" for c in PROD}}

# ---------- 1. VIF screening (iterative to VIF < 5) ----------
# VIF < 5 is the binding requirement. Within each set of strongly correlated predictors the
# more collinear member is dropped iteratively; where two canopy-property classes are themselves
# cross-correlated (e.g. height vs leaf area/heterogeneity), only the least-redundant survives.
vf = df[ALLP].dropna().copy()
for c in ALLP: vf[c] = (vf[c] - vf[c].mean()) / vf[c].std()
keep = ALLP.copy(); vif_log = []
while True:
    vifs = [variance_inflation_factor(vf[keep].values, i) for i in range(len(keep))]
    vmax = max(vifs); jmax = int(np.argmax(vifs))
    if vmax < 5 or len(keep) <= 4: break
    vif_log.append(dict(dropped=keep[jmax], vif=round(vmax, 2), block=blk[keep[jmax]], n_remaining=len(keep)-1))
    keep.pop(jmax)
vif_final = pd.DataFrame({"predictor": keep, "block": [blk[c] for c in keep],
    "VIF": [round(variance_inflation_factor(vf[keep].values, i), 2) for i in range(len(keep))]})
vif_final.to_csv(os.path.join(OUT, "vif_screening.csv"), index=False)
pd.DataFrame(vif_log).to_csv(os.path.join(OUT, "vif_dropped.csv"), index=False)
S_k=[c for c in keep if blk[c]=="structure"]; P_k=[c for c in keep if blk[c]=="spectral"]
D_k=[c for c in keep if blk[c]=="dynamics"];  G_k=[c for c in keep if blk[c]=="productivity"]
mx = max(variance_inflation_factor(vf[keep].values, i) for i in range(len(keep)))
print("VIF-retained (<5):")
for nm,l in [("structure",S_k),("spectral",P_k),("dynamics",D_k),("productivity",G_k)]: print(f"  {nm:12s}: {l}")
print("  max VIF in final set:", round(mx, 2))

# ---------- helpers (OLS + domain fixed) ----------
def fit(resp, preds, d, cluster=False):
    f = f"{resp} ~ C(domain)" + "".join(" + " + p for p in preds)
    if cluster:
        return smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": d["siteID"]})
    return smf.ols(f, d).fit()
def lrt(full, red, dfk):
    chi2 = 2.0 * (full.llf - red.llf)
    return chi2, dfk, (st.chi2.sf(chi2, dfk) if dfk > 0 and chi2 > 0 else np.nan)

# ---------- 2+3. nested models, LRT, variance partitioning ----------
need = list(dict.fromkeys(["domain","siteID","lon","lat"] + keep))
nested_rows, vp_rows, moran_rows, coeff_rows = [], [], [], []
LOGRESP = {"Hill_q1", "Hill_q2"}   # log-transform (satisfies homoscedasticity/normality); LCBD kept raw
for rp, lab in RESP.items():
    d = df[[rp] + need].dropna().copy()
    for c in keep: d[c] = (d[c] - d[c].mean()) / d[c].std()
    d = d.rename(columns={rp: "y"})
    if rp in LOGRESP: d["y"] = np.log(d["y"])
    m0 = fit("y", [], d)                              # domain baseline
    m1 = fit("y", S_k + P_k, d)                       # State
    m2 = fit("y", S_k + P_k + D_k, d)                 # + Dynamics
    m3 = fit("y", S_k + P_k + D_k + G_k, d)           # + Function
    R2dom = m0.rsquared
    for tag, res in [("M0_domain", m0), ("M1_state", m1), ("M2_+dynamics", m2), ("M3_+function", m3)]:
        nested_rows.append(dict(response=lab, model=tag, n=len(d), k=int(res.df_model),
            R2=res.rsquared, R2_adj=res.rsquared_adj, R2_beyond_domain=res.rsquared-R2dom,
            AIC=res.aic, BIC=res.bic, llf=res.llf))
    c12 = lrt(m2, m1, len(D_k)); c23 = lrt(m3, m2, len(G_k))
    nested_rows.append(dict(response=lab, model="LRT_M1->M2(dynamics)", n=len(d),
        dR2=m2.rsquared-m1.rsquared, dAIC=m2.aic-m1.aic, chi2=c12[0], df=c12[1], p=c12[2]))
    nested_rows.append(dict(response=lab, model="LRT_M2->M3(function)", n=len(d),
        dR2=m3.rsquared-m2.rsquared, dAIC=m3.aic-m2.aic, chi2=c23[0], df=c23[1], p=c23[2]))
    # variance partition (unique = R2_full - R2_full_without_block; all include domain)
    blocks = {"structure":S_k, "spectral":P_k, "dynamics":D_k, "productivity":G_k}
    R2full = m3.rsquared; uniq = {}
    for bn, bl in blocks.items():
        rest = [p for p in (S_k+P_k+D_k+G_k) if p not in bl]
        uniq[bn] = R2full - fit("y", rest, d).rsquared
    vp_rows.append(dict(response=lab, R2_full=R2full, R2_domain=R2dom,
        R2_RS=R2full-R2dom, **{f"unique_{k}":v for k,v in uniq.items()},
        shared_RS=(R2full-R2dom)-sum(uniq.values())))
    # cluster-robust coefficients for full M3
    m3c = fit("y", S_k+P_k+D_k+G_k, d, cluster=True)
    for p in (S_k+P_k+D_k+G_k):
        coeff_rows.append(dict(response=lab, predictor=p, block=blk[p],
            beta=m3c.params[p], se_cluster=m3c.bse[p], p_cluster=m3c.pvalues[p]))
    # ---------- 5. Moran's I on M3 residuals ----------
    resid = m3.resid.values
    coords = d[["lon","lat"]].values; k = 8
    nn = NearestNeighbors(n_neighbors=k+1).fit(coords); _, nbr = nn.kneighbors(coords)
    nbr = nbr[:, 1:]; n = len(resid); z = resid - resid.mean()
    moran = lambda zv: (zv * zv[nbr].mean(axis=1)).sum() / (zv**2).sum()
    I_obs = moran(z)
    perm = np.array([moran(np.random.permutation(z)) for _ in range(999)])
    p_perm = (np.sum(perm >= I_obs) + 1) / 1000.0
    # spatial robustness: add smooth spatial trend (lon,lat polynomial), re-partition
    d2 = d.copy(); d2["lon2"]=d2.lon**2; d2["lat2"]=d2.lat**2; d2["lonlat"]=d2.lon*d2.lat
    SP = ["lon","lat","lon2","lat2","lonlat"]
    m3s = fit("y", S_k+P_k+D_k+G_k+SP, d2); m2s = fit("y", S_k+P_k+SP+D_k, d2)  # ref for dyn increment
    uniq_dyn_sp = m3s.rsquared - fit("y", [p for p in S_k+P_k+D_k+G_k if p not in D_k]+SP, d2).rsquared
    uniq_prod_sp = m3s.rsquared - fit("y", S_k+P_k+D_k+SP, d2).rsquared
    moran_rows.append(dict(response=lab, morans_I=I_obs, exp_I=-1/(n-1), p_perm=p_perm, n=n, k=k,
        unique_dynamics_spatial=uniq_dyn_sp, unique_productivity_spatial=uniq_prod_sp))

pd.DataFrame(nested_rows).to_csv(os.path.join(OUT, "nested_models.csv"), index=False)
pd.DataFrame(vp_rows).to_csv(os.path.join(OUT, "variance_partition.csv"), index=False)
pd.DataFrame(moran_rows).to_csv(os.path.join(OUT, "morans_I.csv"), index=False)
pd.DataFrame(coeff_rows).to_csv(os.path.join(OUT, "coeffs_m3_clustered.csv"), index=False)

# ---------- 4. FDR on O4 interactions (DHI excluded) ----------
inter = []
sr = pd.read_csv(os.path.join(R, "obj2_severity_recency.csv"))
inter += [dict(family="disturbance_moderation", response=r.response, term=r.term, p=r.p)
          for r in sr[sr.test.str.contains("x_RS|inter|recency", case=False, na=False)].itertuples()]
am = pd.read_csv(os.path.join(R, "stand_age_moderation_gami.csv"))
inter += [dict(family="standage_moderation", response=r.response, term=r.rs, p=r.p_int) for r in am.itertuples()]
he = pd.read_csv(os.path.join(R, "obj2_heterogeneity.csv"))
if "test" in he:
    inter += [dict(family="landuse_moderation", response=r.response, term=r.term, p=r.p)
              for r in he[he.test.str.contains("x|mod|inter", case=False, na=False)].itertuples()]
fdr = pd.DataFrame(inter).dropna(subset=["p"])
fdr = fdr[~fdr.term.astype(str).str.contains("dhi", case=False, na=False)].reset_index(drop=True)
if len(fdr):
    rej, padj, _, _ = multipletests(fdr.p.values, alpha=0.05, method="fdr_bh")
    fdr["p_fdr"] = padj; fdr["sig_fdr"] = rej
fdr.to_csv(os.path.join(OUT, "interaction_fdr.csv"), index=False)

# ---------- console ----------
pd.set_option("display.width", 200)
print("\n=== NESTED (OLS + domain; R2 / AIC / LRT) ===")
print(pd.DataFrame(nested_rows).round(4).to_string(index=False))
print("\n=== VARIANCE PARTITION (unique/shared beyond domain) ===")
print(pd.DataFrame(vp_rows).round(3).to_string(index=False))
print("\n=== MORAN'S I (M3 residuals) + spatial-robust block uniques ===")
print(pd.DataFrame(moran_rows).round(4).to_string(index=False))
print(f"\n=== INTERACTION FDR: {int(fdr.sig_fdr.sum()) if len(fdr) else 0}/{len(fdr)} survive BH ===")
