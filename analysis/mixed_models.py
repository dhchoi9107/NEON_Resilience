"""
Mixed-Effects Models: RS Diversity → Taxonomic Diversity
=========================================================
Linear mixed models with site nested in domain as random effects.

Models:
  1. Univariate: response ~ predictor + (1|domain/site)
  2. Multivariate: response ~ top predictors + (1|domain/site)

Responses: richness, lcbd_bray
Predictors: structural (13 metrics), spectral (10 metrics)

All predictors z-score standardized within the pooled dataset.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import zscore
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

from site_config import SITE_DOMAIN

OUT_DIR = Path("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs")
RESULT_DIR = Path("E:/neon_lidar/model_results")

DOMAIN_NAMES = {
    "D01": "Northeast", "D02": "Mid-Atlantic", "D03": "Southeast",
    "D05": "Great Lakes", "D07": "Appalachians", "D08": "Ozarks",
    "D10": "Rockies", "D16": "Pacific NW", "D17": "Pacific SW",
}

# ── Load & merge data ─────────────────────────────────────────────────────

print("Loading data...")

# Taxonomic
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness"]]

lcbd_src = pd.read_csv("E:/neon_lidar/model_results/plot_level_dbh10.csv",
                        usecols=["siteID", "plotID", "lcbd_bray"])
lcbd_src = lcbd_src.drop_duplicates("plotID")

taxo = alpha.merge(lcbd_src, on=["siteID", "plotID"], how="inner")
taxo["domain"] = taxo["siteID"].map(SITE_DOMAIN)

# Structural
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

fsd_taxo = fsd.merge(taxo[["siteID", "plotID", "richness", "lcbd_bray"]],
                      on=["siteID", "plotID"], how="inner", suffixes=("_fsd", ""))
print(f"  Structural × Taxo: {len(fsd_taxo)} rows, "
      f"{fsd_taxo['siteID'].nunique()} sites, "
      f"{fsd_taxo.groupby(['siteID','fsd_year']).ngroups} site-years")

# Spectral
vi = pd.read_csv("E:/neon_lidar/spectral_diversity/plot_spectral_1m.csv")
vi = vi[vi["grain_m"] == 1].copy()
vi["domain"] = vi["siteID"].map(SITE_DOMAIN)

spec_cols = [
    "NDVI_mean", "EVI_mean", "ARVI_mean", "PRI_mean", "SAVI_mean",
    "LAI_mean", "fPAR_mean",
    "NDVI_sd", "EVI_sd", "PRI_sd",
]
spec_labels = {
    "NDVI_mean": "NDVI", "EVI_mean": "EVI", "ARVI_mean": "ARVI",
    "PRI_mean": "PRI", "SAVI_mean": "SAVI",
    "LAI_mean": "LAI_opt", "fPAR_mean": "fPAR",
    "NDVI_sd": "NDVI_SD", "EVI_sd": "EVI_SD", "PRI_sd": "PRI_SD",
}

vi_taxo = vi.merge(taxo[["siteID", "plotID", "richness", "lcbd_bray"]],
                    on=["siteID", "plotID"], how="inner")
print(f"  Spectral × Taxo: {len(vi_taxo)} rows, "
      f"{vi_taxo['siteID'].nunique()} sites, "
      f"{vi_taxo.groupby(['siteID','year']).ngroups} site-years")

# ── Z-score standardization ───────────────────────────────────────────────

print("\nStandardizing predictors (z-score)...")

for col in struct_cols:
    if col in fsd_taxo.columns:
        vals = fsd_taxo[col].dropna()
        if vals.std() > 1e-10:
            fsd_taxo[f"z_{struct_labels[col]}"] = (
                (fsd_taxo[col] - vals.mean()) / vals.std()
            )

for col in spec_cols:
    if col in vi_taxo.columns:
        vals = vi_taxo[col].dropna()
        if vals.std() > 1e-10:
            vi_taxo[f"z_{spec_labels[col]}"] = (
                (vi_taxo[col] - vals.mean()) / vals.std()
            )

# Also standardize responses for comparability
for df, resp_cols in [(fsd_taxo, ["richness", "lcbd_bray"]),
                       (vi_taxo, ["richness", "lcbd_bray"])]:
    for col in resp_cols:
        vals = df[col].dropna()
        if vals.std() > 1e-10:
            df[f"z_{col}"] = (df[col] - vals.mean()) / vals.std()

print("  Done.")

# ── Mixed model engine ────────────────────────────────────────────────────

def run_mixed_model(df, response, predictor, site_col="siteID", domain_col="domain"):
    """Run LMM: response ~ predictor + C(domain) + (1|site).

    Domain as fixed effect (9 levels, too few for random).
    Site as random intercept.

    R² marginal includes predictor + domain fixed effects.
    R² predictor isolates the predictor contribution only.
    """
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

    # Extract predictor coefficient
    beta = result.fe_params[predictor]
    se = result.bse_fe[predictor]
    tval = result.tvalues[predictor]
    pval = result.pvalues[predictor]

    # Variance components
    try:
        cov = result.cov_re
        if hasattr(cov, 'iloc') and cov.size > 0:
            var_site = float(cov.iloc[0, 0])
        elif hasattr(cov, 'size') and cov.size > 0:
            var_site = float(cov)
        else:
            var_site = 0.0
    except Exception:
        var_site = 0.0
    var_site = max(var_site, 0.0)
    var_resid = result.scale

    # R² marginal: variance explained by ALL fixed effects (predictor + domain)
    # Use design matrix directly to compute fitted values
    X_design = np.asarray(result.model.exog)
    fe = np.asarray(result.fe_params)
    y_hat_all = X_design @ fe
    var_f_all = np.var(y_hat_all)

    # R² predictor only: variance explained by predictor alone
    # (partial effect, controlling for domain)
    var_f_pred = beta ** 2 * np.var(sub[predictor].values)

    total_var = var_f_all + var_site + var_resid
    r2_marginal = var_f_all / total_var if total_var > 0 else 0
    r2_conditional = (var_f_all + var_site) / total_var if total_var > 0 else 0
    r2_predictor = var_f_pred / total_var if total_var > 0 else 0

    # Domain effect: how much does adding domain improve over predictor alone?
    var_f_domain = var_f_all - var_f_pred
    r2_domain = var_f_domain / total_var if total_var > 0 else 0

    icc_site = var_site / (var_site + var_resid) if (var_site + var_resid) > 0 else 0

    return {
        "response": response,
        "predictor": predictor,
        "n": len(sub),
        "n_domains": sub[domain_col].nunique(),
        "n_sites": sub[site_col].nunique(),
        "beta": beta,
        "se": se,
        "t": tval,
        "p": pval,
        "aic": result.aic,
        "bic": result.bic,
        "var_site": var_site,
        "var_resid": var_resid,
        "icc_site": icc_site,
        "r2_predictor": r2_predictor,
        "r2_domain": r2_domain,
        "r2_marginal": r2_marginal,
        "r2_conditional": r2_conditional,
        "converged": result.converged,
    }


def run_multivariate_mixed(df, response, predictors, site_col="siteID", domain_col="domain"):
    """Run LMM: response ~ pred1 + pred2 + ... + C(domain) + (1|site)."""
    cols = [response] + predictors + [site_col, domain_col]
    sub = df[cols].dropna()
    if len(sub) < 20 or sub[site_col].nunique() < 3 or sub[domain_col].nunique() < 2:
        return None

    formula = f"{response} ~ " + " + ".join(predictors) + f" + C({domain_col})"
    try:
        model = smf.mixedlm(formula, data=sub, groups=sub[site_col])
        result = model.fit(reml=True, method="lbfgs")
    except Exception:
        return None

    if not result.converged:
        return None

    # Extract predictor coefficients
    coefs = {}
    for pred in predictors:
        coefs[pred] = {
            "beta": result.fe_params[pred],
            "se": result.bse_fe[pred],
            "t": result.tvalues[pred],
            "p": result.pvalues[pred],
        }

    try:
        cov = result.cov_re
        if hasattr(cov, 'iloc') and cov.size > 0:
            var_site = float(cov.iloc[0, 0])
        elif hasattr(cov, 'size') and cov.size > 0:
            var_site = float(cov)
        else:
            var_site = 0.0
    except Exception:
        var_site = 0.0
    var_site = max(var_site, 0.0)
    var_resid = result.scale

    # R² marginal (all fixed effects)
    X_design = np.asarray(result.model.exog)
    fe = np.asarray(result.fe_params)
    y_hat_all = X_design @ fe
    var_f_all = np.var(y_hat_all)

    # R² from predictors only (partial)
    X = sub[predictors].values
    betas = np.array([result.fe_params[p] for p in predictors])
    var_f_pred = np.var(X @ betas)

    total = var_f_all + var_site + var_resid
    r2_m = var_f_all / total if total > 0 else 0
    r2_c = (var_f_all + var_site) / total if total > 0 else 0
    r2_pred = var_f_pred / total if total > 0 else 0
    r2_domain = (var_f_all - var_f_pred) / total if total > 0 else 0

    icc_site = var_site / (var_site + var_resid) if (var_site + var_resid) > 0 else 0

    return {
        "response": response,
        "predictors": predictors,
        "n": len(sub),
        "n_domains": sub[domain_col].nunique(),
        "n_sites": sub[site_col].nunique(),
        "coefs": coefs,
        "aic": result.aic,
        "bic": result.bic,
        "var_site": var_site,
        "var_resid": var_resid,
        "icc_site": icc_site,
        "r2_predictor": r2_pred,
        "r2_domain": r2_domain,
        "r2_marginal": r2_m,
        "r2_conditional": r2_c,
    }


# ── Run univariate models ─────────────────────────────────────────────────

print("\n" + "="*80)
print("UNIVARIATE MIXED MODELS: response ~ predictor + C(domain) + (1|site)")
print("="*80)

responses = ["z_richness", "z_lcbd_bray"]
resp_display = {"z_richness": "Richness", "z_lcbd_bray": "LCBD"}

all_univariate = []

# Structural
print("\n--- Structural predictors ---")
struct_z_cols = [f"z_{struct_labels[c]}" for c in struct_cols if f"z_{struct_labels[c]}" in fsd_taxo.columns]

for resp in responses:
    for zcol in struct_z_cols:
        res = run_mixed_model(fsd_taxo, resp, zcol)
        if res:
            res["category"] = "Structural"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_display[resp]
            all_univariate.append(res)

# Spectral
print("--- Spectral predictors ---")
spec_z_cols = [f"z_{spec_labels[c]}" for c in spec_cols if f"z_{spec_labels[c]}" in vi_taxo.columns]

for resp in responses:
    for zcol in spec_z_cols:
        res = run_mixed_model(vi_taxo, resp, zcol)
        if res:
            res["category"] = "Spectral"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_display[resp]
            all_univariate.append(res)

uni_df = pd.DataFrame(all_univariate)

# Print summary table
print(f"\n{'Resp':<10} {'Category':<12} {'Predictor':<12} {'beta':>7} {'SE':>7} "
      f"{'t':>7} {'p':>9} {'R2pred':>7} {'R2dom':>6} {'R2m':>6} {'R2c':>6} {'ICC_s':>6}")
print("-" * 115)

for _, row in uni_df.sort_values(["resp_label", "category", "p"]).iterrows():
    sig = "***" if row["p"] < 0.001 else "** " if row["p"] < 0.01 else "*  " if row["p"] < 0.05 else "   "
    print(f"{row['resp_label']:<10} {row['category']:<12} {row['pred_label']:<12} "
          f"{row['beta']:+7.3f} {row['se']:7.3f} {row['t']:+7.2f} {row['p']:9.2e} "
          f"{row['r2_predictor']:7.3f} {row['r2_domain']:6.3f} "
          f"{row['r2_marginal']:6.3f} {row['r2_conditional']:6.3f} "
          f"{row['icc_site']:6.3f} {sig}")

# Save
uni_df.to_csv(RESULT_DIR / "mixed_model_univariate.csv", index=False)
print(f"\nSaved: mixed_model_univariate.csv ({len(uni_df)} models)")

# ── Run multivariate models ───────────────────────────────────────────────

print("\n" + "="*80)
print("MULTIVARIATE MIXED MODELS")
print("="*80)

# Select top predictors per category (based on univariate p-values)
multi_results = []

for resp in responses:
    resp_label = resp_display[resp]

    # Top structural (by p-value)
    sub = uni_df[(uni_df["response"] == resp) & (uni_df["category"] == "Structural")]
    top_struct = sub.nsmallest(5, "p")["predictor"].tolist()

    if len(top_struct) >= 2:
        print(f"\n{resp_label} ~ Structural (top 5): {[c.replace('z_','') for c in top_struct]}")
        res = run_multivariate_mixed(fsd_taxo, resp, top_struct)
        if res:
            res["category"] = "Structural"
            res["resp_label"] = resp_label
            multi_results.append(res)
            print(f"  R²pred={res['r2_predictor']:.3f}, R²dom={res['r2_domain']:.3f}, "
                  f"R²m={res['r2_marginal']:.3f}, R²c={res['r2_conditional']:.3f}, "
                  f"ICC_site={res['icc_site']:.3f}")
            for pred, c in res["coefs"].items():
                sig = "***" if c["p"] < 0.001 else "** " if c["p"] < 0.01 else "*  " if c["p"] < 0.05 else "   "
                print(f"    {pred.replace('z_',''):<12} beta={c['beta']:+.3f} "
                      f"(SE={c['se']:.3f}, t={c['t']:+.2f}, p={c['p']:.2e}) {sig}")

    # Top spectral
    sub = uni_df[(uni_df["response"] == resp) & (uni_df["category"] == "Spectral")]
    top_spec = sub.nsmallest(5, "p")["predictor"].tolist()

    if len(top_spec) >= 2:
        print(f"\n{resp_label} ~ Spectral (top 5): {[c.replace('z_','') for c in top_spec]}")
        res = run_multivariate_mixed(vi_taxo, resp, top_spec)
        if res:
            res["category"] = "Spectral"
            res["resp_label"] = resp_label
            multi_results.append(res)
            print(f"  R²pred={res['r2_predictor']:.3f}, R²dom={res['r2_domain']:.3f}, "
                  f"R²m={res['r2_marginal']:.3f}, R²c={res['r2_conditional']:.3f}, "
                  f"ICC_site={res['icc_site']:.3f}")
            for pred, c in res["coefs"].items():
                sig = "***" if c["p"] < 0.001 else "** " if c["p"] < 0.01 else "*  " if c["p"] < 0.05 else "   "
                print(f"    {pred.replace('z_',''):<12} beta={c['beta']:+.3f} "
                      f"(SE={c['se']:.3f}, t={c['t']:+.2f}, p={c['p']:.2e}) {sig}")

    # Combined structural + spectral (merge datasets)
    # Need overlapping plots
    merged = fsd_taxo.merge(
        vi_taxo[["siteID", "plotID", "year"] +
                [c for c in vi_taxo.columns if c.startswith("z_")]],
        left_on=["siteID", "plotID", "fsd_year"],
        right_on=["siteID", "plotID", "year"],
        how="inner", suffixes=("", "_spec")
    )

    if len(merged) > 50:
        # Top 3 structural + top 3 spectral
        top_s = [c for c in top_struct[:3] if c in merged.columns]
        top_v = [c for c in top_spec[:3] if c in merged.columns and c not in top_s]
        combined_preds = top_s + top_v

        if len(combined_preds) >= 3:
            print(f"\n{resp_label} ~ Combined (struct+spec): {[c.replace('z_','') for c in combined_preds]}")
            res = run_multivariate_mixed(merged, resp, combined_preds)
            if res:
                res["category"] = "Combined"
                res["resp_label"] = resp_label
                multi_results.append(res)
                print(f"  n={res['n']}, R²pred={res['r2_predictor']:.3f}, R²dom={res['r2_domain']:.3f}, "
                      f"R²m={res['r2_marginal']:.3f}, R²c={res['r2_conditional']:.3f}, "
                      f"ICC_site={res['icc_site']:.3f}")
                for pred, c in res["coefs"].items():
                    sig = "***" if c["p"] < 0.001 else "** " if c["p"] < 0.01 else "*  " if c["p"] < 0.05 else "   "
                    print(f"    {pred.replace('z_',''):<12} beta={c['beta']:+.3f} "
                          f"(SE={c['se']:.3f}, t={c['t']:+.2f}, p={c['p']:.2e}) {sig}")

# ── Visualization ──────────────────────────────────────────────────────────

print("\n" + "="*80)
print("GENERATING FIGURES")
print("="*80)

# Figure 1: Forest plot of univariate coefficients
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
fig.suptitle("Mixed Model Coefficients (Univariate)\n"
             "response ~ z(predictor) + (1|domain/site)\n"
             "Standardized coefficients with 95% CI",
             fontsize=13, fontweight="bold")

panels = [
    ("Richness", "Structural"),
    ("Richness", "Spectral"),
    ("LCBD", "Structural"),
    ("LCBD", "Spectral"),
]

for idx, (resp, cat) in enumerate(panels):
    ax = axes.flatten()[idx]
    sub = uni_df[(uni_df["resp_label"] == resp) & (uni_df["category"] == cat)].copy()
    sub = sub.sort_values("beta")

    y_pos = range(len(sub))
    colors = ["#d73027" if p < 0.001 else "#fc8d59" if p < 0.01
              else "#fee090" if p < 0.05 else "#cccccc"
              for p in sub["p"]]

    ax.barh(list(y_pos), sub["beta"], xerr=1.96 * sub["se"],
            color=colors, edgecolor="black", linewidth=0.5,
            capsize=3, height=0.7)
    ax.axvline(x=0, color="black", linewidth=0.8, linestyle="-")
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(sub["pred_label"], fontsize=9)
    ax.set_xlabel("Standardized coefficient (beta)")
    ax.set_title(f"{resp} ~ {cat}\n"
                 f"R²pred range: {sub['r2_predictor'].min():.3f}-{sub['r2_predictor'].max():.3f}",
                 fontsize=10, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    # Add R²pred annotation
    for i, (_, row) in enumerate(sub.iterrows()):
        ax.text(ax.get_xlim()[1] * 0.95, i,
                f"R²p={row['r2_predictor']:.3f}",
                va="center", ha="right", fontsize=7, style="italic")

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#d73027", label="p < 0.001"),
    Patch(facecolor="#fc8d59", label="p < 0.01"),
    Patch(facecolor="#fee090", label="p < 0.05"),
    Patch(facecolor="#cccccc", label="p >= 0.05"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=4,
           fontsize=10, bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.savefig(OUT_DIR / "mixed_model_forest_plot.png", dpi=150, bbox_inches="tight")
print("Saved: mixed_model_forest_plot.png")
plt.close()

# Figure 2: R² comparison bar chart
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
fig2.suptitle("Variance Explained: Marginal vs Conditional R²\n"
              "(Marginal = fixed effects only, Conditional = fixed + random)",
              fontsize=12, fontweight="bold")

for ax_idx, resp in enumerate(["Richness", "LCBD"]):
    ax = axes2[ax_idx]
    sub = uni_df[uni_df["resp_label"] == resp].copy()
    sub = sub.sort_values("r2_marginal", ascending=True)

    y = range(len(sub))
    labels = [f"{r['pred_label']} ({r['category'][:3]})" for _, r in sub.iterrows()]

    ax.barh(list(y), sub["r2_conditional"], color="#4393c3", alpha=0.5,
            label="R² conditional", height=0.7)
    ax.barh(list(y), sub["r2_marginal"], color="#d6604d",
            label="R² marginal", height=0.7)

    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("R²")
    ax.set_title(resp, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="x", alpha=0.3)

plt.tight_layout()
fig2.savefig(OUT_DIR / "mixed_model_r2_comparison.png", dpi=150, bbox_inches="tight")
print("Saved: mixed_model_r2_comparison.png")
plt.close()

# Figure 3: Variance partitioning — predictor vs domain vs site vs residual
fig3, axes3 = plt.subplots(1, 2, figsize=(16, 8))
fig3.suptitle("Variance Partitioning: Predictor + Domain (fixed) + Site (random) + Residual\n"
              "response ~ z(predictor) + C(domain) + (1|site)",
              fontsize=12, fontweight="bold")

for ax_idx, resp in enumerate(["Richness", "LCBD"]):
    ax = axes3[ax_idx]
    vp = uni_df[uni_df["resp_label"] == resp].copy()
    vp = vp.sort_values("r2_predictor", ascending=True)

    y = range(len(vp))
    labels = [f"{r['pred_label']} ({r['category'][:3]})" for _, r in vp.iterrows()]

    # Stacked bars: predictor + domain + site + residual = 1.0
    r2_resid = 1.0 - vp["r2_marginal"] - vp["icc_site"] * (1 - vp["r2_marginal"])

    ax.barh(list(y), vp["r2_predictor"], color="#d73027",
            edgecolor="black", linewidth=0.3, height=0.7, label="Predictor")
    ax.barh(list(y), vp["r2_domain"], left=vp["r2_predictor"],
            color="#fc8d59", edgecolor="black", linewidth=0.3, height=0.7,
            label="Domain (fixed)")

    # Site random effect share
    site_share = vp["r2_conditional"] - vp["r2_marginal"]
    ax.barh(list(y), site_share, left=vp["r2_marginal"],
            color="#4393c3", edgecolor="black", linewidth=0.3, height=0.7,
            label="Site (random)")

    # Residual
    resid_share = 1.0 - vp["r2_conditional"]
    ax.barh(list(y), resid_share, left=vp["r2_conditional"],
            color="#cccccc", edgecolor="black", linewidth=0.3, height=0.7,
            label="Residual")

    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Proportion of variance")
    ax.set_title(resp, fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    ax.set_xlim(0, 1.0)

plt.tight_layout()
fig3.savefig(OUT_DIR / "mixed_model_variance_partition.png", dpi=150, bbox_inches="tight")
print("Saved: mixed_model_variance_partition.png")
plt.close()

# ── Summary ────────────────────────────────────────────────────────────────

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

for resp in ["Richness", "LCBD"]:
    sub = uni_df[uni_df["resp_label"] == resp].sort_values("p")
    print(f"\n{resp} — Top 5 predictors (by p-value):")
    for i, (_, r) in enumerate(sub.head(5).iterrows()):
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        print(f"  {i+1}. {r['pred_label']:12s} ({r['category']:10s}): "
              f"beta={r['beta']:+.3f}, R²pred={r['r2_predictor']:.3f}, "
              f"R²dom={r['r2_domain']:.3f}, R²c={r['r2_conditional']:.3f}, p={r['p']:.2e} {sig}")

    sig_count = (sub["p"] < 0.05).sum()
    print(f"  → {sig_count}/{len(sub)} predictors significant (p<0.05)")
    print(f"  → Mean R²_predictor = {sub['r2_predictor'].mean():.3f} "
          f"(predictor alone)")
    print(f"  → Mean R²_domain   = {sub['r2_domain'].mean():.3f} "
          f"(domain fixed effect)")
    print(f"  → Mean ICC_site    = {sub['icc_site'].mean():.3f} "
          f"({sub['icc_site'].mean()*100:.0f}% between sites)")

# ── Interaction models: predictor × domain ─────────────────────────────────

print("\n" + "="*80)
print("INTERACTION MODELS: response ~ predictor * C(domain) + (1|site)")
print("="*80)

DOMAINS_ORDERED = ["D01", "D02", "D03", "D05", "D07", "D08", "D10", "D16", "D17"]

def run_interaction_model(df, response, predictor, site_col="siteID", domain_col="domain"):
    """Run LMM with interaction: response ~ predictor * C(domain) + (1|site).

    Returns main effect, interaction significance, and per-domain slopes.
    """
    required = [response, predictor, site_col, domain_col]
    sub = df[required].dropna()
    if len(sub) < 30 or sub[site_col].nunique() < 3 or sub[domain_col].nunique() < 3:
        return None

    # Main effect model (no interaction)
    formula_main = f"{response} ~ {predictor} + C({domain_col})"
    # Interaction model
    formula_int = f"{response} ~ {predictor} * C({domain_col})"

    try:
        model_main = smf.mixedlm(formula_main, data=sub, groups=sub[site_col])
        model_int = smf.mixedlm(formula_int, data=sub, groups=sub[site_col])

        # Try lbfgs first, fallback to powell
        for method in ["lbfgs", "powell", "cg"]:
            try:
                result_main = model_main.fit(reml=False, method=method)
                result_int = model_int.fit(reml=False, method=method)
                if result_int.converged:
                    break
            except Exception:
                continue
        else:
            return None
    except Exception:
        return None

    if not result_int.converged:
        return None

    # LRT: interaction model vs main effect model
    ll_main = result_main.llf
    ll_int = result_int.llf
    n_int_params = len(result_int.fe_params) - len(result_main.fe_params)
    from scipy.stats import chi2
    lrt_stat = 2 * (ll_int - ll_main)
    lrt_p = chi2.sf(lrt_stat, df=n_int_params) if n_int_params > 0 else 1.0

    # R² comparison
    X_main = np.asarray(result_main.model.exog)
    fe_main = np.asarray(result_main.fe_params)
    var_main = np.var(X_main @ fe_main)

    X_int = np.asarray(result_int.model.exog)
    fe_int = np.asarray(result_int.fe_params)
    var_int = np.var(X_int @ fe_int)

    try:
        cov = result_int.cov_re
        var_site = float(cov.iloc[0, 0]) if hasattr(cov, 'iloc') and cov.size > 0 else 0.0
    except:
        var_site = 0.0
    var_site = max(var_site, 0.0)
    var_resid = result_int.scale

    total = var_int + var_site + var_resid
    r2_int = var_int / total if total > 0 else 0
    try:
        cov_main = result_main.cov_re
        var_site_main = float(cov_main.iloc[0,0]) if hasattr(cov_main,'iloc') and cov_main.size > 0 else 0.0
    except Exception:
        var_site_main = 0.0
    var_site_main = max(var_site_main, 0.0)
    total_main = var_main + var_site_main + result_main.scale
    r2_main = var_main / total_main if total_main > 0 else 0
    r2_interaction_gain = r2_int - r2_main

    # Per-domain slopes (reference domain + interaction terms)
    ref_domain = sorted(sub[domain_col].unique())[0]
    base_slope = result_int.fe_params[predictor]

    domain_slopes = {}
    for d in sorted(sub[domain_col].unique()):
        if d == ref_domain:
            domain_slopes[d] = base_slope
        else:
            int_key = f"{predictor}:C({domain_col})[T.{d}]"
            if int_key in result_int.fe_params:
                domain_slopes[d] = base_slope + result_int.fe_params[int_key]
            else:
                domain_slopes[d] = np.nan

    return {
        "response": response,
        "predictor": predictor,
        "n": len(sub),
        "lrt_stat": lrt_stat,
        "lrt_p": lrt_p,
        "n_int_params": n_int_params,
        "r2_main": r2_main,
        "r2_interaction": r2_int,
        "r2_gain": r2_interaction_gain,
        "aic_main": result_main.aic,
        "aic_int": result_int.aic,
        "domain_slopes": domain_slopes,
    }


# Run interaction models for top predictors
top_struct = ["z_Deep_Gap", "z_Rumple", "z_Ht_Ratio", "z_FHD", "z_Vert_CV",
              "z_LAI", "z_Gap_Frac", "z_VCI"]
top_spec = ["z_NDVI", "z_EVI", "z_ARVI", "z_PRI", "z_SAVI", "z_fPAR"]

interaction_results = []

for resp in responses:
    resp_label = resp_display[resp]

    for zcol in top_struct:
        if zcol not in fsd_taxo.columns:
            continue
        res = run_interaction_model(fsd_taxo, resp, zcol)
        if res:
            res["category"] = "Structural"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_label
            interaction_results.append(res)

    for zcol in top_spec:
        if zcol not in vi_taxo.columns:
            continue
        res = run_interaction_model(vi_taxo, resp, zcol)
        if res:
            res["category"] = "Spectral"
            res["pred_label"] = zcol.replace("z_", "")
            res["resp_label"] = resp_label
            interaction_results.append(res)

# Print interaction results
print(f"\n{'Resp':<10} {'Cat':<6} {'Predictor':<12} {'LRT_p':>10} {'R2_main':>8} "
      f"{'R2_int':>7} {'R2_gain':>8} {'AIC_m':>8} {'AIC_int':>8} {'Sig?'}")
print("-" * 100)

for r in sorted(interaction_results, key=lambda x: (x["resp_label"], x["lrt_p"])):
    sig = "***" if r["lrt_p"] < 0.001 else "** " if r["lrt_p"] < 0.01 else "*  " if r["lrt_p"] < 0.05 else "   "
    print(f"{r['resp_label']:<10} {r['category'][:5]:<6} {r['pred_label']:<12} "
          f"{r['lrt_p']:10.2e} {r['r2_main']:8.3f} {r['r2_interaction']:7.3f} "
          f"{r['r2_gain']:+8.3f} {r['aic_main']:8.1f} {r['aic_int']:8.1f} {sig}")

# Domain-specific slopes heatmap
print("\n--- Per-domain slopes (predictor × domain interaction) ---")
print(f"\n{'Resp':<10} {'Predictor':<12}", end="")
for d in DOMAINS_ORDERED:
    print(f" {d:>6}", end="")
print(f"  {'LRT_p':>8}")

for r in sorted(interaction_results, key=lambda x: (x["resp_label"], x["category"], x["pred_label"])):
    print(f"{r['resp_label']:<10} {r['pred_label']:<12}", end="")
    for d in DOMAINS_ORDERED:
        slope = r["domain_slopes"].get(d, np.nan)
        if np.isnan(slope):
            print(f"    -- ", end="")
        else:
            print(f" {slope:+.3f}", end="")
    sig = "***" if r["lrt_p"] < 0.001 else "**" if r["lrt_p"] < 0.01 else "*" if r["lrt_p"] < 0.05 else ""
    print(f"  {r['lrt_p']:8.2e} {sig}")

# Heatmap figure: domain-specific slopes
print("\n--- Generating interaction heatmap ---")

for resp_label in ["Richness", "LCBD"]:
    int_sub = [r for r in interaction_results if r["resp_label"] == resp_label]
    if not int_sub:
        continue

    pred_labels = [r["pred_label"] for r in int_sub]
    n_pred = len(pred_labels)
    n_dom = len(DOMAINS_ORDERED)

    slope_mat = np.full((n_dom, n_pred), np.nan)
    for j, r in enumerate(int_sub):
        for i, d in enumerate(DOMAINS_ORDERED):
            slope_mat[i, j] = r["domain_slopes"].get(d, np.nan)

    fig_h = max(4, 0.5 * n_dom + 2)
    fig_w = max(8, 0.8 * n_pred + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = np.nanmax(np.abs(slope_mat)) * 1.0
    im = ax.imshow(slope_mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    for i in range(n_dom):
        for j in range(n_pred):
            if np.isnan(slope_mat[i, j]):
                continue
            sv = slope_mat[i, j]
            color = "white" if abs(sv) > vmax * 0.6 else "black"
            ax.text(j, i, f"{sv:+.2f}", ha="center", va="center",
                    fontsize=8, color=color, fontweight="bold" if abs(sv) > vmax * 0.4 else "normal")

    domain_labels = [f"{d} {DOMAIN_NAMES.get(d, '')}" for d in DOMAINS_ORDERED]
    ax.set_yticks(range(n_dom))
    ax.set_yticklabels(domain_labels, fontsize=9)
    ax.set_xticks(range(n_pred))

    # Add significance stars to x-labels
    xlabels = []
    for r in int_sub:
        sig = "***" if r["lrt_p"] < 0.001 else "**" if r["lrt_p"] < 0.01 else "*" if r["lrt_p"] < 0.05 else ""
        cat_short = "S" if r["category"] == "Structural" else "V"
        xlabels.append(f"{r['pred_label']}({cat_short}){sig}")
    ax.set_xticklabels(xlabels, fontsize=8, rotation=45, ha="right")

    ax.set_title(f"{resp_label} ~ Predictor × Domain Interaction\n"
                 f"Cell = domain-specific slope (from interaction model)",
                 fontsize=11, fontweight="bold", pad=10)

    plt.colorbar(im, ax=ax, label="Domain-specific slope", shrink=0.8, pad=0.02)
    plt.tight_layout()

    fname = f"mixed_model_interaction_{resp_label.lower()}.png"
    fig.savefig(OUT_DIR / fname, dpi=150, bbox_inches="tight")
    print(f"Saved: {fname}")
    plt.close()

# Save interaction results
int_rows = []
for r in interaction_results:
    row = {k: v for k, v in r.items() if k != "domain_slopes"}
    for d, slope in r["domain_slopes"].items():
        row[f"slope_{d}"] = slope
    int_rows.append(row)
pd.DataFrame(int_rows).to_csv(RESULT_DIR / "mixed_model_interaction.csv", index=False)
print(f"Saved: mixed_model_interaction.csv")

print("\nDone.")
