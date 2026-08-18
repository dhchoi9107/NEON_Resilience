# -*- coding: utf-8 -*-
"""
Model diagnostics for the plot-level framework (the checks the manuscript asserts
in Methods 2.5 but that were never re-derived): predictor collinearity (VIF),
residual heteroscedasticity (Breusch-Pagan), residual normality (Shapiro-Wilk,
plus skew/kurtosis), and residual spatial autocorrelation (Moran's I, k = 8
nearest neighbours, permutation test). Raw vs log scale is compared for the Hill
numbers to document the transformation decision.
Out: results/O0_framework/model_diagnostics.csv
"""
import io, sys, os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8",errors="replace")
import statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D, RES = os.path.join(BASE,"data"), os.path.join(BASE,"papers","UNIFIED","results","O0_framework")
S  = ["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]; P = ["EVI_mean"]
Dn = ["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]
keep = S+P+Dn
RESP = {"Hill_q1":"Hill q1","Hill_q2":"Hill q2","LCBD_turnover_rare":"Turnover","LCBD_nestedness_rare":"Nestedness"}
LOG  = {"Hill_q1","Hill_q2"}

f = pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv")); f = f[f.sample_coverage>=0.9].copy()
ll = pd.read_csv(os.path.join(D,"plot_lonlat_26.csv"))[["plotID","lon","lat"]]
f = f.merge(ll, on="plotID", how="left")

# ---------- 1. VIF on the retained predictor set ----------
print("="*70); print("1. COLLINEARITY (VIF) - manuscript claims all < 5, max 4.0"); print("="*70)
dv = f[keep].dropna().copy()
for c in keep: dv[c] = (dv[c]-dv[c].mean())/dv[c].std()
X = np.column_stack([np.ones(len(dv))] + [dv[c].values for c in keep])
vifs = {c: variance_inflation_factor(X, i+1) for i, c in enumerate(keep)}
for c, v in sorted(vifs.items(), key=lambda kv: -kv[1]): print(f"   {c:20s} VIF = {v:.2f}")
print(f"   -> max VIF = {max(vifs.values()):.2f}  (claim: 4.0, all < 5)")

# ---------- 2-4. per-response diagnostics ----------
def moran_I(resid, xy, k=8, nperm=999, seed=11):
    from scipy.spatial import cKDTree
    n = len(resid); tree = cKDTree(xy)
    _, idx = tree.query(xy, k=k+1)          # first neighbour is self
    W = np.zeros((n, n))
    for i in range(n): W[i, idx[i,1:]] = 1.0
    W = W / W.sum(1, keepdims=True)
    z = resid - resid.mean()
    num = z @ (W @ z); den = (z**2).sum()
    I = n/W.sum() * num/den
    rng = np.random.RandomState(seed); cnt = 0
    for _ in range(nperm):
        zp = rng.permutation(z)
        Ip = n/W.sum() * (zp @ (W @ zp))/((zp**2).sum())
        if abs(Ip) >= abs(I): cnt += 1
    return I, (cnt+1)/(nperm+1)

rows = []
print("\n" + "="*70); print("2. RESIDUAL DIAGNOSTICS (full model: domain + structure + spectral + dynamics)"); print("="*70)
for rr, lab in RESP.items():
    d = f[[rr,"domain","siteID","lon","lat"]+keep].dropna().copy()
    for c in keep: d[c] = (d[c]-d[c].mean())/d[c].std()
    for scale in (["raw","log"] if rr in LOG else ["raw"]):
        y = np.log(d[rr]) if scale == "log" else d[rr].astype(float)
        d["y"] = (y-y.mean())/y.std()
        m = smf.ols("y ~ C(domain) + " + " + ".join(keep), d).fit()
        r = m.resid.values
        bp = het_breuschpagan(r, m.model.exog)[1]
        sh = stats.shapiro(r if len(r) <= 5000 else r[:5000])[1]
        xy = d[["lon","lat"]].values
        I, pI = moran_I(r, xy)
        rows.append(dict(response=lab, scale=scale, n=len(d),
                         breusch_pagan_p=round(bp,4), shapiro_p=round(sh,4),
                         skew=round(float(stats.skew(r)),3), kurtosis=round(float(stats.kurtosis(r)),3),
                         morans_I=round(I,4), morans_p=round(pI,4)))
        print(f"   {lab:11s} [{scale:3s}] n={len(d)}  BP p={bp:.3f}  Shapiro p={sh:.3f}  "
              f"skew={stats.skew(r):+.2f}  Moran I={I:+.3f} (p={pI:.3f})")

out = pd.DataFrame(rows)
out.insert(0, "max_VIF", round(max(vifs.values()),2))
out.to_csv(os.path.join(RES,"model_diagnostics.csv"), index=False)
print("\nsaved -> model_diagnostics.csv")
