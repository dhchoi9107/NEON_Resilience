# -*- coding: utf-8 -*-
"""
Code-result consistency audit.
Recomputes the manuscript's quantitative claims directly from the source data and
compares them with (a) the stored result CSVs and (b) the numbers written in the
manuscript text. Prints PASS/FAIL per claim. Nothing is taken on trust.
"""
import os, re, warnings, numpy as np, pandas as pd, scipy.stats as st
warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
RES  = os.path.join(BASE, "papers", "UNIFIED", "results", "O0_framework")
D    = os.path.join(BASE, "data")
MS   = os.path.join(BASE, "papers", "UNIFIED", "MANUSCRIPT_v4_final.md")
TXT  = open(MS, encoding="utf-8").read()

S  = ["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]
P  = ["EVI_mean"]
Dn = ["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]
keep = S + P + Dn
RESP = {"Hill_q1":"Hill q1","Hill_q2":"Hill q2",
        "LCBD_turnover_rare":"LCBD turnover","LCBD_nestedness_rare":"LCBD nestedness"}
LOG = {"Hill_q1","Hill_q2"}

f = pd.read_csv(os.path.join(D,"FINAL_v2_pooled_26.csv"))
f = f[f.sample_coverage >= 0.9].copy()

n_pass = n_fail = 0
def chk(label, got, want, tol=0.0015, unit=""):
    global n_pass, n_fail
    ok = abs(got-want) <= tol
    n_pass, n_fail = n_pass+ok, n_fail+(not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: recomputed={got:.4f}{unit} manuscript={want}{unit}")
    return ok
def chk_int(label, got, want):
    global n_pass, n_fail
    ok = got == want
    n_pass, n_fail = n_pass+ok, n_fail+(not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: recomputed={got} manuscript={want}")
def in_text(label, s):
    global n_pass, n_fail
    ok = s in TXT
    n_pass, n_fail = n_pass+ok, n_fail+(not ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] text contains: {label}")

def prep(resp, cols):
    d = f[[resp,"domain","siteID"]+cols].dropna().copy()
    y = np.log(d[resp]) if resp in LOG else d[resp].astype(float)
    d["y"] = (y-y.mean())/y.std()
    for c in cols: d[c] = (d[c]-d[c].mean())/d[c].std()
    return d
def r2(d, pr):
    fml = "y ~ C(domain)" + ("" if not pr else " + " + " + ".join(pr))
    return smf.ols(fml, d).fit().rsquared

print("="*72); print("A. SAMPLE SIZES"); print("="*72)
d_a = prep("Hill_q1", keep); d_b = prep("LCBD_turnover_rare", keep)
chk_int("alpha plots (n)", len(d_a), 644)
chk_int("beta plots (n)",  len(d_b), 579)
chk_int("sites",  f.siteID.nunique(), 26)
chk_int("domains", f.domain.nunique(), 12)

print("="*72); print("B. VARIANCE PARTITION / NESTED R2  (recomputed from raw data)"); print("="*72)
vp = pd.read_csv(os.path.join(RES,"variance_partition_sd.csv")).set_index("response")
nm = pd.read_csv(os.path.join(RES,"nested_models_sd.csv")).set_index("response")
want_dom = {"Hill_q1":0.49,"Hill_q2":0.45,"LCBD_turnover_rare":0.58,"LCBD_nestedness_rare":0.33}
want_seq = {"Hill_q1":0.042,"Hill_q2":0.039,"LCBD_turnover_rare":0.013,"LCBD_nestedness_rare":0.056}
want_uni_str = {"Hill_q1":0.081,"Hill_q2":0.073}
for rr, lab in RESP.items():
    d = prep(rr, keep)
    dom  = r2(d, [])
    m1   = r2(d, S+P) - dom          # state beyond domain
    m2   = r2(d, S+P+Dn) - dom       # + dynamics beyond domain
    seq  = m2 - m1                   # sequential dynamics increment
    uS   = m2 - (r2(d, P+Dn)  - dom)
    uP   = m2 - (r2(d, S+Dn)  - dom)
    uDn  = m2 - (r2(d, S+P)   - dom)
    print(f"- {lab}")
    chk("  R2 domain", dom, want_dom[rr], tol=0.006)
    chk("  seq dynamics increment", seq, want_seq[rr], tol=0.0015)
    chk("  CSV agreement: unique_dynamics", uDn, float(vp.loc[lab,"unique_dynamics"]), tol=0.0015)
    chk("  CSV agreement: unique_structure", uS, float(vp.loc[lab,"unique_structure"]), tol=0.0015)
    chk("  CSV agreement: unique_spectral", uP, float(vp.loc[lab,"unique_spectral"]), tol=0.0015)
    chk("  CSV agreement: R2_M2_beyond", m2, float(nm.loc[lab,"R2_M2_beyond"]), tol=0.0015)
    if rr in want_uni_str: chk("  unique structure (text)", uS, want_uni_str[rr], tol=0.0015)
    if rr == "LCBD_nestedness_rare": chk("  unique dynamics (text 0.056)", uDn, 0.056, tol=0.0015)

print("="*72); print("C. BAYESIAN COEFFICIENTS (CSV vs manuscript text)"); print("="*72)
bc = pd.read_csv(os.path.join(RES,"bayes_multilevel_coeffs.csv"))
def bget(resp, pred, col="beta"):
    r = bc[(bc.response==resp)&(bc.predictor==pred)]
    return float(r[col].iloc[0]) if len(r) else np.nan
for pred, want in [("Ht_Ratio_trend",-0.31),("Vert_CV_trend",0.20),("LAI_trend",-0.19)]:
    chk(f"nestedness {pred}", bget("LCBD_nestedness_rare",pred), want, tol=0.006)
chk("Hill q1 LAI_mean", bget("Hill_q1","LAI_mean"), 0.41, tol=0.006)
chk("Hill q1 VCI_mean", bget("Hill_q1","VCI_mean"), 0.27, tol=0.006)
chk("Hill q2 LAI_mean", bget("Hill_q2","LAI_mean"), 0.40, tol=0.006)
chk("Hill q2 VCI_mean", bget("Hill_q2","VCI_mean"), 0.26, tol=0.006)
chk("nestedness Vert_CV_mean", bget("LCBD_nestedness_rare","Vert_CV_mean"), -0.12, tol=0.006)
chk("nestedness Rugosity_mean", bget("LCBD_nestedness_rare","Rugosity_mean"), 0.14, tol=0.006)
# turnover: two small credible dynamics trends ~ +0.09
tv = bc[(bc.response=="LCBD_turnover_rare")&(bc.credible)]
chk_int("turnover credible count", len(tv), 2)
print("     turnover credible:", [(r.predictor, round(r.beta,3)) for r in tv.itertuples()])
# alpha dynamics credible (Results 3.4 claims rumple +0.09, ht-ratio +0.11)
chk("Hill q1 Rumple_trend", bget("Hill_q1","Rumple_trend"), 0.09, tol=0.006)
chk("Hill q1 Ht_Ratio_trend", bget("Hill_q1","Ht_Ratio_trend"), 0.11, tol=0.007)

print("="*72); print("D. SITE-LEVEL SPECIES-ENERGY (recomputed)"); print("="*72)
g3 = pd.read_csv(os.path.join(D,"plot_pml_gpp_ts_26.csv"))[["plotID","pml_gpp"]]
gm = pd.read_csv(os.path.join(D,"plot_modis_gpp_26.csv"))[["plotID","modis_gpp"]]
site = (f.merge(g3,on="plotID").merge(gm,on="plotID",how="left")
          .groupby("siteID").agg(Hill_q1=("Hill_q1","mean"),Hill_q2=("Hill_q2","mean"),
                                 pml=("pml_gpp","mean"),modis=("modis_gpp","mean")).reset_index())
tw = pd.read_csv(os.path.join(D,"site_tower_gpp_26.csv"))[["siteID","tower_gpp"]]
site = site.merge(tw,on="siteID",how="left"); tv2 = site.dropna(subset=["tower_gpp"])
chk_int("tower sites", len(tv2), 22)
r,p = st.pearsonr(tv2.tower_gpp, tv2.Hill_q1); chk("tower vs Hill q1 r", r, 0.59, tol=0.006); chk("  p", p, 0.004, tol=0.0006)
r2_,p2_ = st.pearsonr(tv2.tower_gpp, tv2.Hill_q2); chk("tower vs Hill q2 r", r2_, 0.57, tol=0.006); chk("  p", p2_, 0.006, tol=0.0009)
z = (tv2.tower_gpp-tv2.tower_gpp.mean())/tv2.tower_gpp.std()
qp = smf.ols("Hill_q1 ~ z + I(z**2)", tv2.assign(z=z)).fit().pvalues["I(z ** 2)"]
chk("tower quadratic p (no hump)", qp, 0.87, tol=0.01)
rp,pp = st.pearsonr(site.pml, site.Hill_q1);  chk("PML vs Hill q1 r", rp, 0.50, tol=0.006); chk("  p", pp, 0.009, tol=0.001)
sm_ = site.dropna(subset=["modis"]); rm,pm = st.pearsonr(sm_.modis, sm_.Hill_q1)
chk("MODIS cross-check r (Table S9)", rm, 0.60, tol=0.006)
chk("PML vs tower r",   st.pearsonr(tv2.pml, tv2.tower_gpp)[0], 0.86, tol=0.006)
chk("MODIS vs tower r (Table S9)", st.pearsonr(tv2.modis, tv2.tower_gpp)[0], 0.85, tol=0.007)

print("="*72); print("E. ICC of satellite GPP (claim: ~0.92)"); print("="*72)
gp = f.merge(g3, on="plotID")
grand = gp.pml_gpp.mean(); sw = gp.groupby("siteID").pml_gpp
ss_between = (sw.mean().sub(grand).pow(2) * sw.size()).sum()
ss_total   = ((gp.pml_gpp-grand)**2).sum()
chk("between-site share of plot GPP variance", ss_between/ss_total, 0.92, tol=0.02)

print("="*72); print("F. INTERACTION SCREEN (O4)"); print("="*72)
o4 = pd.read_csv(os.path.join(RES,"o4_interactions_pml.csv"))
chk_int("screen terms", len(o4), 144)
chk_int("FDR-surviving", int((o4.q<0.05).sum()), 5)
fam = o4[o4.q<0.05].family.value_counts().to_dict()
chk_int("  disturbance", fam.get("disturbance",0), 2)
chk_int("  land use",    fam.get("land use",0), 3)
ss_ = pd.read_csv(os.path.join(RES,"simple_slopes_pml.csv"))
gA = ss_[(ss_.rp=="Hill q1")&(ss_.rs=="pml_gpp")&(ss_.ctx=="stand_age_gami")].sort_values("ctx_z")
gB = ss_[(ss_.rp=="Hill q1")&(ss_.rs=="VCI_mean")&(ss_.ctx=="severity")].sort_values("ctx_z")
for v,w in zip(gA.slope,[0.04,0.12,0.20]): chk("PMLxage slope", v, w, tol=0.006)
for v,w in zip(gB.slope,[0.39,0.30,0.21]): chk("VCIxseverity slope", v, w, tol=0.006)

print("="*72); print("G. TEXT-ONLY CLAIMS (presence check)"); print("="*72)
for lab,s in [("Mantel r +0.32/+0.30","Spearman r = +0.32 and +0.30"),
              ("8,813 pairs","8,813 within-site pairs"),
              ("pseudo-R2 0.11","pseudo-R² = 0.11"),
              ("BACI -1.05 p=0.069","β = −1.05, p = 0.069"),
              ("stand age n=639","n = 639 aged plots"),
              ("press -0.029 / pulse ~0.000","press median −0.029, 23 sites")]:
    in_text(lab, s)

print("="*72)
print(f"SUMMARY: {n_pass} PASS / {n_fail} FAIL")
print("="*72)
