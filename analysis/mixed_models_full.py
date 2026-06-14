"""
Full Mixed-Effects Models: All Components → Taxonomic Diversity & Productivity
===============================================================================
Integrates structural, spectral, functional diversity, environmental
heterogeneity, and productivity (DHI) into comprehensive mixed models.

Two analysis levels:
  A) Plot-level: structural + spectral + functional → richness / LCBD
  B) Site-year level: all components → beta diversity, gamma richness

Model form: response ~ z(predictors) + C(domain) + (1|site)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from site_config import SITE_DOMAIN

OUT_DIR = Path("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs")
RESULT_DIR = Path("E:/neon_lidar/model_results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DOMAIN_NAMES = {
    "D01": "Northeast", "D02": "Mid-Atlantic", "D03": "Southeast",
    "D05": "Great Lakes", "D07": "Appalachians", "D08": "Ozarks",
    "D10": "Rockies", "D16": "Pacific NW", "D17": "Pacific SW",
}
DOMAINS_ORDERED = ["D01", "D02", "D03", "D05", "D07", "D08", "D10", "D16", "D17"]


# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

print("Loading data...")

# ── Taxonomic diversity is already in plot_level_complete.csv ──
# (richness, shannon, simpson, lcbd_bray columns from previous pipeline)

# ── Structural (plot-level) ──
fsd = pd.read_csv("E:/neon_lidar/model_results/plot_level_complete.csv")
fsd["domain"] = fsd["siteID"].map(SITE_DOMAIN)

struct_cols = [
    "FHD_mean", "rumple_mean", "mean_max_canopy_ht_mean", "max_canopy_ht_mean",
    "vert_sd_mean", "vertCV_mean", "LAI_mean", "GC_mean",
    "deepgap_fraction_mean", "VCI_mean", "top_rugosity_mean",
    "q95_mean", "HeightRatio_mean",
]
struct_labels = {
    "FHD_mean": "FHD", "rumple_mean": "Rumple",
    "mean_max_canopy_ht_mean": "Canopy_Ht", "max_canopy_ht_mean": "Max_Ht",
    "vert_sd_mean": "Vert_SD", "vertCV_mean": "Vert_CV",
    "LAI_mean": "LAI", "GC_mean": "Gap_Frac",
    "deepgap_fraction_mean": "Deep_Gap", "VCI_mean": "VCI",
    "top_rugosity_mean": "Rugosity", "q95_mean": "Q95",
    "HeightRatio_mean": "Ht_Ratio",
}

# ── Spectral (plot-level, 1m) ──
vi = pd.read_csv("E:/neon_lidar/spectral_diversity/plot_spectral_1m.csv")
vi = vi[vi["grain_m"] == 1].copy()
vi["domain"] = vi["siteID"].map(SITE_DOMAIN)

spec_cols = [
    "NDVI_mean", "EVI_mean", "ARVI_mean", "PRI_mean", "SAVI_mean",
    "LAI_mean", "fPAR_mean", "NDVI_sd", "EVI_sd", "PRI_sd",
]
spec_labels = {
    "NDVI_mean": "NDVI", "EVI_mean": "EVI", "ARVI_mean": "ARVI",
    "PRI_mean": "PRI", "SAVI_mean": "SAVI",
    "LAI_mean": "LAI_opt", "fPAR_mean": "fPAR",
    "NDVI_sd": "NDVI_SD", "EVI_sd": "EVI_SD", "PRI_sd": "PRI_SD",
}

# ── Functional diversity (plot-level) ──
func = pd.read_csv(RESULT_DIR / "plot_functional_diversity.csv")
func["domain"] = func["siteID"].map(SITE_DOMAIN)

func_cols = ["func_FRic", "func_FDiv", "func_FEve", "func_RaoQ"]
func_labels = {
    "func_FRic": "FRic", "func_FDiv": "FDiv",
    "func_FEve": "FEve", "func_RaoQ": "Func_RaoQ",
}

# ── Site-year level assembled data ──
assembled = pd.read_csv(RESULT_DIR / "assembled_data.csv")
assembled["domain"] = assembled["siteID"].map(SITE_DOMAIN)

site_year_preds = {
    # Spectral diversity
    "rao_q": "Spec_RaoQ", "spectral_cv": "Spec_CV",
    "spectral_shannon": "Spec_Shannon",
    "spectral_FRic": "Spec_FRic", "spectral_FDiv": "Spec_FDiv",
    # Env heterogeneity
    "chm_cv": "CHM_CV", "chm_mean": "CHM_Mean", "chm_iqr": "CHM_IQR",
    "Sa_100m": "Sa_100m", "Sa_500m": "Sa_500m", "Sa_1000m": "Sa_1000m",
    "Sq_500m": "Sq_500m", "Ssk_500m": "Ssk_500m", "Sku_500m": "Sku_500m",
    # Functional diversity
    "func_FRic": "Func_FRic_sy", "func_FDiv": "Func_FDiv_sy",
    "func_RaoQ": "Func_RaoQ_sy",
    # Productivity
    "cumulative_mean": "DHI_Cum", "minimum_mean": "DHI_Min",
    "variation_mean": "DHI_Var", "trend_mean": "Ht_Trend",
}
site_year_responses = {
    "alpha_richness_mean": "Alpha_Rich",
    "alpha_shannon_mean": "Alpha_Shannon",
    "bray_mean": "Beta_Bray",
    "gamma_richness": "Gamma_Rich",
}

print(f"  Structural: {len(fsd)} plot-years")
print(f"  Spectral: {len(vi)} plot-years")
print(f"  Functional: {len(func)} plot-years")
print(f"  Assembled (site-year): {len(assembled)} rows")


# ═══════════════════════════════════════════════════════════════════════════
# MERGE & STANDARDIZE
# ═══════════════════════════════════════════════════════════════════════════

print("\nMerging plot-level data...")

# A) Structural — fsd already contains richness/shannon/simpson/lcbd_bray
fsd_taxo = fsd[fsd["richness"].notna()].copy()
print(f"  Structural (with taxo): {len(fsd_taxo)} rows, "
      f"{fsd_taxo['siteID'].nunique()} sites")

# B) Structural + Functional (year-matched)
fsd_func_taxo = fsd.merge(func, left_on=["siteID", "plotID", "fsd_year"],
                           right_on=["siteID", "plotID", "year"],
                           how="inner", suffixes=("", "_f"))
fsd_func_taxo = fsd_func_taxo[fsd_func_taxo["richness"].notna()].copy()
print(f"  Struct × Func (with taxo): {len(fsd_func_taxo)} rows")

# C) Spectral × Taxonomic (need to bring in diversity from fsd)
# Use pooled alpha since spectral has different years
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness", "shannon", "simpson"]]
lcbd_path = RESULT_DIR / "plot_level_dbh10.csv"
if lcbd_path.exists():
    lcbd = pd.read_csv(lcbd_path, usecols=["siteID", "plotID", "lcbd_bray"])
    lcbd = lcbd.drop_duplicates("plotID")
    taxo = alpha.merge(lcbd, on=["siteID", "plotID"], how="inner")
else:
    taxo = alpha.copy()
    taxo["lcbd_bray"] = np.nan

vi_taxo = vi.merge(taxo, on=["siteID", "plotID"], how="inner")
print(f"  Spectral × Taxo: {len(vi_taxo)} rows")

# ── Z-score standardization ──
print("\nStandardizing predictors...")

def z_standardize(df, cols, labels):
    """Z-score standardize columns, adding z_Label columns."""
    for col in cols:
        if col not in df.columns:
            continue
        label = labels[col]
        vals = df[col].dropna()
        if len(vals) > 10 and vals.std() > 1e-10:
            df[f"z_{label}"] = (df[col] - vals.mean()) / vals.std()

# Plot-level
z_standardize(fsd_taxo, struct_cols, struct_labels)
z_standardize(fsd_func_taxo, struct_cols, struct_labels)
z_standardize(fsd_func_taxo, func_cols, func_labels)
z_standardize(vi_taxo, spec_cols, spec_labels)

# Standardize responses
for df in [fsd_taxo, fsd_func_taxo, vi_taxo]:
    for col in ["richness", "shannon", "simpson", "lcbd_bray"]:
        if col in df.columns:
            vals = df[col].dropna()
            if vals.std() > 1e-10:
                df[f"z_{col}"] = (df[col] - vals.mean()) / vals.std()

# Site-year level
for col, label in site_year_preds.items():
    if col in assembled.columns:
        vals = assembled[col].dropna()
        if len(vals) > 5 and vals.std() > 1e-10:
            assembled[f"z_{label}"] = (assembled[col] - vals.mean()) / vals.std()

for col, label in site_year_responses.items():
    if col in assembled.columns:
        vals = assembled[col].dropna()
        if len(vals) > 5 and vals.std() > 1e-10:
            assembled[f"z_{label}"] = (assembled[col] - vals.mean()) / vals.std()

print("  Done.")


# ═══════════════════════════════════════════════════════════════════════════
# MODEL ENGINE
# ═══════════════════════════════════════════════════════════════════════════

def run_lmm(df, response, predictor, site_col="siteID", domain_col="domain"):
    """LMM: response ~ predictor + C(domain) + (1|site)."""
    required = [response, predictor, site_col, domain_col]
    sub = df[required].dropna()
    if len(sub) < 20 or sub[site_col].nunique() < 3 or sub[domain_col].nunique() < 2:
        return None

    formula = f"{response} ~ {predictor} + C({domain_col})"
    try:
        model = smf.mixedlm(formula, data=sub, groups=sub[site_col])
        result = model.fit(reml=True, method="lbfgs")
    except Exception:
        return None

    if not result.converged:
        return None

    beta = result.fe_params[predictor]
    se = result.bse_fe[predictor]
    pval = result.pvalues[predictor]

    # Variance components
    try:
        cov = result.cov_re
        var_site = float(cov.iloc[0, 0]) if hasattr(cov, 'iloc') and cov.size > 0 else float(cov) if hasattr(cov, 'size') and cov.size > 0 else 0.0
    except Exception:
        var_site = 0.0
    var_site = max(var_site, 0.0)
    var_resid = result.scale

    X = np.asarray(result.model.exog)
    fe = np.asarray(result.fe_params)
    y_hat = X @ fe
    var_f_all = np.var(y_hat)
    var_f_pred = beta ** 2 * np.var(sub[predictor].values)

    total = var_f_all + var_site + var_resid
    r2m = var_f_all / total if total > 0 else 0
    r2c = (var_f_all + var_site) / total if total > 0 else 0
    r2p = var_f_pred / total if total > 0 else 0
    r2d = (var_f_all - var_f_pred) / total if total > 0 else 0
    icc = var_site / (var_site + var_resid) if (var_site + var_resid) > 0 else 0

    return {
        "response": response, "predictor": predictor,
        "n": len(sub), "n_sites": sub[site_col].nunique(),
        "n_domains": sub[domain_col].nunique(),
        "beta": beta, "se": se, "t": result.tvalues[predictor], "p": pval,
        "aic": result.aic, "bic": result.bic,
        "var_site": var_site, "var_resid": var_resid, "icc_site": icc,
        "r2_predictor": r2p, "r2_domain": r2d,
        "r2_marginal": r2m, "r2_conditional": r2c,
    }


def run_multi_lmm(df, response, predictors, site_col="siteID", domain_col="domain"):
    """Multivariate LMM: response ~ pred1 + pred2 + ... + C(domain) + (1|site)."""
    cols = [response] + predictors + [site_col, domain_col]
    sub = df[[c for c in cols if c in df.columns]].dropna()
    avail_preds = [p for p in predictors if p in sub.columns]
    if len(sub) < 20 or sub[site_col].nunique() < 3 or len(avail_preds) < 2:
        return None

    formula = f"{response} ~ " + " + ".join(avail_preds) + f" + C({domain_col})"
    for method in ["lbfgs", "powell", "cg"]:
        try:
            model = smf.mixedlm(formula, data=sub, groups=sub[site_col])
            result = model.fit(reml=True, method=method)
            if result.converged:
                break
        except Exception:
            continue
    else:
        return None

    if not result.converged:
        return None

    coefs = {}
    for pred in avail_preds:
        coefs[pred] = {
            "beta": result.fe_params[pred], "se": result.bse_fe[pred],
            "t": result.tvalues[pred], "p": result.pvalues[pred],
        }

    try:
        cov = result.cov_re
        var_site = float(cov.iloc[0, 0]) if hasattr(cov, 'iloc') and cov.size > 0 else 0.0
    except Exception:
        var_site = 0.0
    var_site = max(var_site, 0.0)
    var_resid = result.scale

    X = np.asarray(result.model.exog)
    fe = np.asarray(result.fe_params)
    var_f_all = np.var(X @ fe)

    betas = np.array([result.fe_params[p] for p in avail_preds])
    var_f_pred = np.var(sub[avail_preds].values @ betas)

    total = var_f_all + var_site + var_resid
    r2m = var_f_all / total if total > 0 else 0
    r2c = (var_f_all + var_site) / total if total > 0 else 0
    r2p = var_f_pred / total if total > 0 else 0
    icc = var_site / (var_site + var_resid) if (var_site + var_resid) > 0 else 0

    return {
        "response": response, "predictors": avail_preds, "coefs": coefs,
        "n": len(sub), "n_sites": sub[site_col].nunique(),
        "aic": result.aic, "bic": result.bic,
        "var_site": var_site, "var_resid": var_resid, "icc_site": icc,
        "r2_predictor": r2p, "r2_domain": (var_f_all - var_f_pred) / total if total > 0 else 0,
        "r2_marginal": r2m, "r2_conditional": r2c,
    }


# ═══════════════════════════════════════════════════════════════════════════
# A) PLOT-LEVEL UNIVARIATE MODELS
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("A) PLOT-LEVEL UNIVARIATE MODELS")
print("   response ~ z(predictor) + C(domain) + (1|site)")
print("=" * 80)

responses_plot = {
    "z_richness": "Richness",
    "z_shannon": "Shannon",
    "z_lcbd_bray": "LCBD",
}

all_uni = []

# Structural predictors
print("\n--- Structural predictors ---")
struct_z = [f"z_{struct_labels[c]}" for c in struct_cols
            if f"z_{struct_labels[c]}" in fsd_taxo.columns]

for resp, resp_label in responses_plot.items():
    if resp not in fsd_taxo.columns:
        continue
    for zcol in struct_z:
        res = run_lmm(fsd_taxo, resp, zcol)
        if res:
            res["category"] = "Structural"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_label
            res["level"] = "plot"
            all_uni.append(res)

# Spectral predictors
print("--- Spectral predictors ---")
spec_z = [f"z_{spec_labels[c]}" for c in spec_cols
          if f"z_{spec_labels[c]}" in vi_taxo.columns]

for resp, resp_label in responses_plot.items():
    if resp not in vi_taxo.columns:
        continue
    for zcol in spec_z:
        res = run_lmm(vi_taxo, resp, zcol)
        if res:
            res["category"] = "Spectral"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_label
            res["level"] = "plot"
            all_uni.append(res)

# Functional predictors
print("--- Functional predictors ---")
func_z = [f"z_{func_labels[c]}" for c in func_cols
          if f"z_{func_labels[c]}" in fsd_func_taxo.columns]

for resp, resp_label in responses_plot.items():
    if resp not in fsd_func_taxo.columns:
        continue
    for zcol in func_z:
        res = run_lmm(fsd_func_taxo, resp, zcol)
        if res:
            res["category"] = "Functional"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_label
            res["level"] = "plot"
            all_uni.append(res)

uni_df = pd.DataFrame(all_uni)

# Print summary
print(f"\n{'Resp':<10} {'Cat':<12} {'Predictor':<12} {'beta':>7} {'p':>10} "
      f"{'R2pred':>7} {'R2dom':>6} {'R2m':>6} {'R2c':>6} {'n':>5}")
print("-" * 100)

for _, r in uni_df.sort_values(["resp_label", "category", "p"]).iterrows():
    sig = "***" if r["p"] < 0.001 else "** " if r["p"] < 0.01 else "*  " if r["p"] < 0.05 else "   "
    print(f"{r['resp_label']:<10} {r['category']:<12} {r['pred_label']:<12} "
          f"{r['beta']:+7.3f} {r['p']:10.2e} "
          f"{r['r2_predictor']:7.3f} {r['r2_domain']:6.3f} "
          f"{r['r2_marginal']:6.3f} {r['r2_conditional']:6.3f} "
          f"{r['n']:5d} {sig}")


# ═══════════════════════════════════════════════════════════════════════════
# B) SITE-YEAR LEVEL MODELS
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("B) SITE-YEAR LEVEL MODELS (assembled data)")
print("   response ~ z(predictor) + C(domain) + (1|site)")
print("=" * 80)

sy_uni = []
sy_preds_z = [f"z_{v}" for v in site_year_preds.values()
              if f"z_{v}" in assembled.columns]
sy_resp_z = [f"z_{v}" for v in site_year_responses.values()
             if f"z_{v}" in assembled.columns]

for resp_col, resp_label in site_year_responses.items():
    z_resp = f"z_{resp_label}"
    if z_resp not in assembled.columns:
        continue
    for pred_col, pred_label in site_year_preds.items():
        z_pred = f"z_{pred_label}"
        if z_pred not in assembled.columns:
            continue
        res = run_lmm(assembled, z_resp, z_pred)
        if res:
            res["pred_label"] = pred_label
            res["resp_label"] = resp_label
            res["level"] = "site-year"
            # Categorize
            if pred_label.startswith("Spec_"):
                res["category"] = "Spectral"
            elif pred_label.startswith("Func_"):
                res["category"] = "Functional"
            elif pred_label.startswith("DHI_") or pred_label == "Ht_Trend":
                res["category"] = "Productivity"
            else:
                res["category"] = "Env_Het"
            sy_uni.append(res)

sy_df = pd.DataFrame(sy_uni)

if len(sy_df) > 0:
    print(f"\n{'Resp':<14} {'Cat':<12} {'Predictor':<14} {'beta':>7} {'p':>10} "
          f"{'R2pred':>7} {'R2m':>6} {'R2c':>6} {'n':>4}")
    print("-" * 100)

    for _, r in sy_df.sort_values(["resp_label", "category", "p"]).iterrows():
        sig = "***" if r["p"] < 0.001 else "** " if r["p"] < 0.01 else "*  " if r["p"] < 0.05 else "   "
        print(f"{r['resp_label']:<14} {r['category']:<12} {r['pred_label']:<14} "
              f"{r['beta']:+7.3f} {r['p']:10.2e} "
              f"{r['r2_predictor']:7.3f} {r['r2_marginal']:6.3f} "
              f"{r['r2_conditional']:6.3f} {r['n']:4d} {sig}")


# ═══════════════════════════════════════════════════════════════════════════
# C) MULTIVARIATE MODELS (best predictors per category)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("C) MULTIVARIATE MODELS (top predictors)")
print("=" * 80)

multi_results = []

# C1: Plot-level — Structural + Functional → Richness/Shannon
for resp, resp_label in [("z_richness", "Richness"), ("z_shannon", "Shannon")]:
    if resp not in fsd_func_taxo.columns:
        continue

    # Get top structural + functional predictors from univariate
    sub_s = uni_df[(uni_df["resp_label"] == resp_label) & (uni_df["category"] == "Structural")]
    sub_f = uni_df[(uni_df["resp_label"] == resp_label) & (uni_df["category"] == "Functional")]

    top_s = sub_s.nsmallest(3, "p")["predictor"].tolist() if len(sub_s) > 0 else []
    top_f = sub_f.nsmallest(2, "p")["predictor"].tolist() if len(sub_f) > 0 else []

    combined = [c for c in top_s + top_f if c in fsd_func_taxo.columns]
    if len(combined) >= 3:
        print(f"\n{resp_label} ~ Structural + Functional: {[c.replace('z_','') for c in combined]}")
        res = run_multi_lmm(fsd_func_taxo, resp, combined)
        if res:
            res["resp_label"] = resp_label
            res["model_name"] = "Struct+Func"
            multi_results.append(res)
            print(f"  n={res['n']}, R2p={res['r2_predictor']:.3f}, R2m={res['r2_marginal']:.3f}, "
                  f"R2c={res['r2_conditional']:.3f}")
            for p, c in res["coefs"].items():
                sig = "***" if c["p"] < 0.001 else "**" if c["p"] < 0.01 else "*" if c["p"] < 0.05 else ""
                print(f"    {p.replace('z_',''):<14} beta={c['beta']:+.3f} p={c['p']:.2e} {sig}")

# C2: Site-year — Spectral + EnvHet + Functional → Beta diversity
for resp_col, resp_label in [("z_Beta_Bray", "Beta_Bray"), ("z_Alpha_Rich", "Alpha_Rich"),
                               ("z_Gamma_Rich", "Gamma_Rich")]:
    if resp_col not in assembled.columns:
        continue

    sub_sy = sy_df[sy_df["resp_label"] == resp_label.replace("z_", "")]
    if len(sub_sy) == 0:
        # Try matching
        sub_sy = sy_df[sy_df["resp_label"] == resp_label]
    if len(sub_sy) < 3:
        continue

    top_preds = sub_sy.nsmallest(5, "p")["predictor"].tolist()
    top_preds = [p for p in top_preds if p in assembled.columns]

    if len(top_preds) >= 3:
        print(f"\n{resp_label} ~ Top 5 (site-year): {[c.replace('z_','') for c in top_preds]}")
        res = run_multi_lmm(assembled, resp_col, top_preds)
        if res:
            res["resp_label"] = resp_label
            res["model_name"] = "AllComponents"
            multi_results.append(res)
            print(f"  n={res['n']}, R2p={res['r2_predictor']:.3f}, R2m={res['r2_marginal']:.3f}, "
                  f"R2c={res['r2_conditional']:.3f}")
            for p, c in res["coefs"].items():
                sig = "***" if c["p"] < 0.001 else "**" if c["p"] < 0.01 else "*" if c["p"] < 0.05 else ""
                print(f"    {p.replace('z_',''):<14} beta={c['beta']:+.3f} p={c['p']:.2e} {sig}")

# C3: Productivity model — DHI ~ diversity measures (site-year)
print("\n--- Productivity models: DHI ~ diversity ---")
for dhi_resp in ["z_DHI_Cum", "z_Ht_Trend"]:
    if dhi_resp not in assembled.columns:
        continue
    dhi_label = dhi_resp.replace("z_", "")

    # Try diversity predictors
    div_preds = [p for p in ["z_Spec_RaoQ", "z_Spec_Shannon", "z_CHM_CV",
                              "z_Sa_500m", "z_Func_FRic_sy", "z_Func_RaoQ_sy"]
                 if p in assembled.columns]

    if len(div_preds) >= 3:
        print(f"\n{dhi_label} ~ Diversity predictors: {[c.replace('z_','') for c in div_preds]}")
        res = run_multi_lmm(assembled, dhi_resp, div_preds)
        if res:
            res["resp_label"] = dhi_label
            res["model_name"] = "DHI~Diversity"
            multi_results.append(res)
            print(f"  n={res['n']}, R2p={res['r2_predictor']:.3f}, R2m={res['r2_marginal']:.3f}, "
                  f"R2c={res['r2_conditional']:.3f}")
            for p, c in res["coefs"].items():
                sig = "***" if c["p"] < 0.001 else "**" if c["p"] < 0.01 else "*" if c["p"] < 0.05 else ""
                print(f"    {p.replace('z_',''):<14} beta={c['beta']:+.3f} p={c['p']:.2e} {sig}")


# ═══════════════════════════════════════════════════════════════════════════
# D) FIGURES
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("GENERATING FIGURES")
print("=" * 80)

# Combine univariate results
all_results = pd.concat([uni_df, sy_df], ignore_index=True)

# ── Figure 1: Forest plot — all univariate coefficients by response ──
for resp_label in ["Richness", "Shannon", "LCBD"]:
    sub = all_results[(all_results["resp_label"] == resp_label) &
                       (all_results["level"] == "plot")].copy()
    if sub.empty:
        continue

    sub = sub.sort_values("beta")
    fig, ax = plt.subplots(figsize=(10, max(6, len(sub) * 0.35)))

    y_pos = range(len(sub))
    cat_colors = {"Structural": "#2166ac", "Spectral": "#b2182b", "Functional": "#1b7837"}
    colors = [cat_colors.get(c, "#999999") for c in sub["category"]]

    ax.barh(list(y_pos), sub["beta"], xerr=1.96 * sub["se"],
            color=colors, edgecolor="black", linewidth=0.3,
            capsize=2, height=0.7, alpha=0.8)
    ax.axvline(x=0, color="black", linewidth=0.8)

    labels = [f"{r['pred_label']} ({r['category'][:3]})" for _, r in sub.iterrows()]
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Standardized coefficient (beta ± 95% CI)")
    ax.set_title(f"{resp_label} ~ z(predictor) + C(domain) + (1|site)\n"
                 f"Plot-level models (n = {sub['n'].iloc[0]})",
                 fontsize=11, fontweight="bold")

    # Significance markers
    for i, (_, r) in enumerate(sub.iterrows()):
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        ax.text(ax.get_xlim()[1], i, f" R²p={r['r2_predictor']:.3f} {sig}",
                va="center", fontsize=7)

    from matplotlib.patches import Patch
    legend_elements = [Patch(fc=c, label=l) for l, c in cat_colors.items()]
    ax.legend(handles=legend_elements, fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    fname = f"mixed_full_forest_{resp_label.lower()}.png"
    fig.savefig(OUT_DIR / fname, dpi=150, bbox_inches="tight")
    print(f"Saved: {fname}")
    plt.close()

# ── Figure 2: Site-year level heatmap (predictor × response) ──
if len(sy_df) > 0:
    fig, ax = plt.subplots(figsize=(14, 8))

    resp_labels = sorted(sy_df["resp_label"].unique())
    pred_labels = sorted(sy_df["pred_label"].unique())

    beta_mat = np.full((len(resp_labels), len(pred_labels)), np.nan)
    sig_mat = np.full((len(resp_labels), len(pred_labels)), "", dtype=object)

    for _, r in sy_df.iterrows():
        ri = resp_labels.index(r["resp_label"])
        pi = pred_labels.index(r["pred_label"])
        beta_mat[ri, pi] = r["beta"]
        sig_mat[ri, pi] = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""

    vmax = np.nanmax(np.abs(beta_mat))
    im = ax.imshow(beta_mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    for i in range(len(resp_labels)):
        for j in range(len(pred_labels)):
            if np.isnan(beta_mat[i, j]):
                continue
            color = "white" if abs(beta_mat[i, j]) > vmax * 0.6 else "black"
            ax.text(j, i, f"{beta_mat[i,j]:+.2f}{sig_mat[i,j]}",
                    ha="center", va="center", fontsize=7, color=color)

    ax.set_xticks(range(len(pred_labels)))
    ax.set_xticklabels(pred_labels, fontsize=8, rotation=45, ha="right")
    ax.set_yticks(range(len(resp_labels)))
    ax.set_yticklabels(resp_labels, fontsize=10)
    ax.set_title("Site-Year Level: Standardized Coefficients\n"
                 "response ~ z(predictor) + C(domain) + (1|site)",
                 fontsize=12, fontweight="bold")
    plt.colorbar(im, ax=ax, label="Beta", shrink=0.6)
    plt.tight_layout()

    fig.savefig(OUT_DIR / "mixed_full_siteyear_heatmap.png", dpi=150, bbox_inches="tight")
    print("Saved: mixed_full_siteyear_heatmap.png")
    plt.close()

# ── Figure 3: Variance partitioning (all categories) ──
fig, axes = plt.subplots(1, 3, figsize=(20, 8))
fig.suptitle("Variance Partitioning: Plot-level Models\n"
             "response ~ z(predictor) + C(domain) + (1|site)",
             fontsize=13, fontweight="bold")

for ax_idx, resp in enumerate(["Richness", "Shannon", "LCBD"]):
    ax = axes[ax_idx]
    vp = all_results[(all_results["resp_label"] == resp) &
                      (all_results["level"] == "plot")].copy()
    if vp.empty:
        continue
    vp = vp.sort_values("r2_predictor", ascending=True)

    y = range(len(vp))
    labels = [f"{r['pred_label']} ({r['category'][:3]})" for _, r in vp.iterrows()]

    # Stacked bars
    ax.barh(list(y), vp["r2_predictor"], color="#d73027",
            edgecolor="black", linewidth=0.2, height=0.7, label="Predictor")
    ax.barh(list(y), vp["r2_domain"], left=vp["r2_predictor"],
            color="#fc8d59", edgecolor="black", linewidth=0.2, height=0.7,
            label="Domain")
    site_share = vp["r2_conditional"] - vp["r2_marginal"]
    ax.barh(list(y), site_share, left=vp["r2_marginal"],
            color="#4393c3", edgecolor="black", linewidth=0.2, height=0.7,
            label="Site")
    resid = 1.0 - vp["r2_conditional"]
    ax.barh(list(y), resid, left=vp["r2_conditional"],
            color="#cccccc", edgecolor="black", linewidth=0.2, height=0.7,
            label="Residual")

    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Proportion of variance")
    ax.set_title(resp, fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_xlim(0, 1)
    ax.grid(axis="x", alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT_DIR / "mixed_full_variance_partition.png", dpi=150, bbox_inches="tight")
print("Saved: mixed_full_variance_partition.png")
plt.close()

# ── Figure 4: R² comparison — structural vs spectral vs functional ──
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("R² predictor by Category (Plot-level Models)",
             fontsize=13, fontweight="bold")

for ax_idx, resp in enumerate(["Richness", "Shannon", "LCBD"]):
    ax = axes[ax_idx]
    sub = all_results[(all_results["resp_label"] == resp) &
                       (all_results["level"] == "plot")].copy()
    if sub.empty:
        continue

    cat_colors = {"Structural": "#2166ac", "Spectral": "#b2182b", "Functional": "#1b7837"}
    for cat in ["Structural", "Spectral", "Functional"]:
        cat_sub = sub[sub["category"] == cat].sort_values("r2_predictor", ascending=False)
        if cat_sub.empty:
            continue
        x = range(len(cat_sub))
        bars = ax.bar([xi + list(cat_colors.keys()).index(cat) * 0.25 for xi in x],
                      cat_sub["r2_predictor"], width=0.23,
                      color=cat_colors[cat], label=cat, alpha=0.8)

    ax.set_xlabel("Predictor rank")
    ax.set_ylabel("R² predictor")
    ax.set_title(resp, fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT_DIR / "mixed_full_r2_by_category.png", dpi=150, bbox_inches="tight")
print("Saved: mixed_full_r2_by_category.png")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

all_results.to_csv(RESULT_DIR / "mixed_model_full_univariate.csv", index=False)
print(f"\nSaved: mixed_model_full_univariate.csv ({len(all_results)} models)")

if sy_df is not None and len(sy_df) > 0:
    sy_df.to_csv(RESULT_DIR / "mixed_model_siteyear.csv", index=False)
    print(f"Saved: mixed_model_siteyear.csv ({len(sy_df)} models)")


# ═══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\n--- Plot-level (univariate) ---")
for resp in ["Richness", "Shannon", "LCBD"]:
    sub = uni_df[uni_df["resp_label"] == resp].sort_values("p")
    if sub.empty:
        continue
    sig_n = (sub["p"] < 0.05).sum()
    print(f"\n{resp}: {sig_n}/{len(sub)} significant (p<0.05)")
    print(f"  Top 5:")
    for i, (_, r) in enumerate(sub.head(5).iterrows()):
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        print(f"    {i+1}. {r['pred_label']:12s} ({r['category']:10s}): "
              f"beta={r['beta']:+.3f}, R2p={r['r2_predictor']:.3f}, p={r['p']:.2e} {sig}")

if len(sy_df) > 0:
    print("\n--- Site-year level ---")
    for resp_label in sorted(sy_df["resp_label"].unique()):
        sub = sy_df[sy_df["resp_label"] == resp_label].sort_values("p")
        sig_n = (sub["p"] < 0.05).sum()
        print(f"\n{resp_label}: {sig_n}/{len(sub)} significant (p<0.05)")
        for i, (_, r) in enumerate(sub.head(3).iterrows()):
            sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
            print(f"    {i+1}. {r['pred_label']:14s} ({r['category']:10s}): "
                  f"beta={r['beta']:+.3f}, R2p={r['r2_predictor']:.3f}, p={r['p']:.2e} {sig}")

if multi_results:
    print("\n--- Multivariate models ---")
    for res in multi_results:
        print(f"\n  {res['resp_label']} ~ {res['model_name']}: "
              f"R2p={res['r2_predictor']:.3f}, R2m={res['r2_marginal']:.3f}, "
              f"R2c={res['r2_conditional']:.3f}, n={res['n']}")

print("\nDone.")
