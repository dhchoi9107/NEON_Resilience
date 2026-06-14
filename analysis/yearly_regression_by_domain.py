"""
Within-Site Regression Grouped by NEON Domain
================================================
Are RS-taxonomy relationships consistent across ecoregions?

Domains in our data:
  D01: Northeast (BART, HARV)
  D02: Mid-Atlantic (BLAN, SCBI, SERC)
  D03: Southeast (JERC, OSBS)
  D05: Great Lakes (STEI, UNDE)   [CHEQ, TREE excluded if no data]
  D07: Appalachians (GRSM, MLBS, ORNL)
  D08: Ozarks/South Central (TALL)
  D10: Central Plains/Rockies (RMNP)
  D16: Pacific Northwest (ABBY, WREF)
  D17: Pacific Southwest (SOAP, TEAK)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from site_config import SITE_DOMAIN

OUT_DIR = Path("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs")

DOMAIN_NAMES = {
    "D01": "Northeast",
    "D02": "Mid-Atlantic",
    "D03": "Southeast",
    "D05": "Great Lakes",
    "D07": "Appalachians",
    "D08": "Ozarks/S.Central",
    "D10": "C.Plains/Rockies",
    "D16": "Pacific NW",
    "D17": "Pacific SW",
}

# ── Load data ───────────────────────────────────────────────────────────────

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
fsd_taxo = fsd.merge(taxo[["siteID", "plotID", "richness", "lcbd_bray"]],
                      on=["siteID", "plotID"], how="inner", suffixes=("_fsd", ""))

# Spectral
vi = pd.read_csv("E:/neon_lidar/spectral_diversity/plot_spectral_1m.csv")
vi = vi[vi["grain_m"] == 1].copy()
vi["domain"] = vi["siteID"].map(SITE_DOMAIN)
vi_taxo = vi.merge(taxo[["siteID", "plotID", "richness", "lcbd_bray"]],
                    on=["siteID", "plotID"], how="inner")

print(f"Structural × Taxo: {len(fsd_taxo)} rows")
print(f"Spectral × Taxo: {len(vi_taxo)} rows")

# ── Regression per site-year ────────────────────────────────────────────────

def site_year_regression(df, x_col, y_col, min_n=8):
    results = []
    for (site, year), sub in df.groupby(["siteID", "fsd_year" if "fsd_year" in df.columns else "year"]):
        sub = sub.dropna(subset=[x_col, y_col])
        if len(sub) < min_n:
            continue
        x, y = sub[x_col].values, sub[y_col].values
        if np.std(x) < 1e-10 or np.std(y) < 1e-10:
            continue
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        domain = sub["domain"].iloc[0]
        results.append({
            "siteID": site, "year": int(year), "domain": domain,
            "n": len(sub), "r": r_value, "r2": r_value**2, "p_value": p_value,
        })
    return pd.DataFrame(results)

# Key predictors
struct_pairs = [
    ("FHD_mean", "richness", "FHD → Richness"),
    ("deepgap_fraction_mean", "richness", "Deep Gap → Richness"),
    ("LAI_mean", "richness", "LAI → Richness"),
    ("vert_sd_mean", "richness", "Vert SD → Richness"),
    ("rumple_mean", "richness", "Rumple → Richness"),
    ("vertCV_mean", "richness", "Vert CV → Richness"),
    ("VCI_mean", "richness", "VCI → Richness"),
    ("FHD_mean", "lcbd_bray", "FHD → LCBD"),
    ("vert_sd_mean", "lcbd_bray", "Vert SD → LCBD"),
    ("LAI_mean", "lcbd_bray", "LAI → LCBD"),
    ("mean_max_canopy_ht_mean", "lcbd_bray", "Canopy Ht → LCBD"),
    ("q95_mean", "lcbd_bray", "Q95 → LCBD"),
]

spec_pairs = [
    ("NDVI_mean", "richness", "NDVI → Richness"),
    ("EVI_sd", "richness", "EVI SD → Richness"),
    ("ARVI_mean", "richness", "ARVI → Richness"),
    ("PRI_mean", "richness", "PRI → Richness"),
    ("NDVI_mean", "lcbd_bray", "NDVI → LCBD"),
    ("EVI_sd", "lcbd_bray", "EVI SD → LCBD"),
    ("ARVI_mean", "lcbd_bray", "ARVI → LCBD"),
    ("LAI_mean", "lcbd_bray", "LAI(opt) → LCBD"),
]

# Run all
all_results = {}

for x_col, y_col, label in struct_pairs:
    if x_col not in fsd_taxo.columns:
        continue
    res = site_year_regression(fsd_taxo, x_col, y_col)
    if len(res) > 0:
        all_results[label] = res

for x_col, y_col, label in spec_pairs:
    df_use = vi_taxo.copy()
    if "fsd_year" not in df_use.columns:
        df_use["fsd_year"] = df_use["year"]
    if x_col not in df_use.columns:
        continue
    res = site_year_regression(df_use, x_col, y_col)
    if len(res) > 0:
        all_results[label] = res

# ── Figure: Domain-level summary ────────────────────────────────────────────

# For each predictor, group results by domain and compute sig rate + mean r
domains_ordered = ["D01", "D02", "D03", "D05", "D07", "D08", "D10", "D16", "D17"]

# Select top predictors for each category
top_predictors = {
    "Richness ~ Structural": ["FHD → Richness", "Vert CV → Richness", "LAI → Richness",
                               "Rumple → Richness", "VCI → Richness"],
    "Richness ~ Spectral": ["NDVI → Richness", "EVI SD → Richness", "ARVI → Richness",
                             "PRI → Richness"],
    "LCBD ~ Structural": ["Vert SD → LCBD", "LAI → LCBD", "FHD → LCBD",
                           "Canopy Ht → LCBD", "Q95 → LCBD"],
    "LCBD ~ Spectral": ["NDVI → LCBD", "EVI SD → LCBD", "ARVI → LCBD",
                         "LAI(opt) → LCBD"],
}

fig, axes = plt.subplots(2, 2, figsize=(20, 16))
fig.suptitle("Within-Site Regression by NEON Domain\n"
             "(bar height = sig rate; color intensity = mean |r|)",
             fontsize=14, fontweight="bold")

cat_colors = {
    "Richness ~ Structural": plt.cm.Blues,
    "Richness ~ Spectral": plt.cm.Oranges,
    "LCBD ~ Structural": plt.cm.Greens,
    "LCBD ~ Spectral": plt.cm.Purples,
}

for ax_idx, (cat, pred_list) in enumerate(top_predictors.items()):
    ax = axes.flatten()[ax_idx]
    cmap = cat_colors[cat]

    # Build grouped data
    group_data = []
    for pred in pred_list:
        if pred not in all_results:
            continue
        res = all_results[pred]
        for domain in domains_ordered:
            d_res = res[res["domain"] == domain]
            if len(d_res) == 0:
                continue
            sig_rate = (d_res["p_value"] < 0.05).mean() * 100
            mean_r = d_res["r"].mean()
            n_sy = len(d_res)
            group_data.append({
                "predictor": pred.split("→")[0].strip(),
                "domain": domain,
                "domain_name": DOMAIN_NAMES.get(domain, domain),
                "sig_rate": sig_rate,
                "mean_r": mean_r,
                "n_site_years": n_sy,
            })

    gdf = pd.DataFrame(group_data)
    if gdf.empty:
        ax.set_visible(False)
        continue

    predictors_in = sorted(gdf["predictor"].unique())
    domains_in = [d for d in domains_ordered if d in gdf["domain"].values]
    n_pred = len(predictors_in)
    n_dom = len(domains_in)

    bar_width = 0.8 / n_pred
    x = np.arange(n_dom)

    for i, pred in enumerate(predictors_in):
        rates = []
        mean_rs = []
        for d in domains_in:
            row = gdf[(gdf["predictor"] == pred) & (gdf["domain"] == d)]
            if len(row) > 0:
                rates.append(row["sig_rate"].values[0])
                mean_rs.append(row["mean_r"].values[0])
            else:
                rates.append(0)
                mean_rs.append(0)

        offset = (i - n_pred / 2 + 0.5) * bar_width
        color_intensity = [min(abs(mr) / 0.4, 1.0) * 0.7 + 0.3 for mr in mean_rs]
        bar_colors = [cmap(ci) for ci in color_intensity]

        bars = ax.bar(x + offset, rates, bar_width * 0.9, color=bar_colors,
                       edgecolor="black", linewidth=0.5, label=pred)

        # Add mean r text
        for j, (rate, mr) in enumerate(zip(rates, mean_rs)):
            if rate > 0:
                ax.text(x[j] + offset, rate + 0.5, f"{mr:+.2f}",
                        ha="center", va="bottom", fontsize=5.5, rotation=90)

    domain_labels = [f"{d}\n{DOMAIN_NAMES.get(d, '')}" for d in domains_in]
    ax.set_xticks(x)
    ax.set_xticklabels(domain_labels, fontsize=8)
    ax.set_ylabel("Significant site-years (%)")
    ax.set_title(cat, fontsize=12, fontweight="bold")
    ax.axhline(y=5, color="red", linestyle="--", alpha=0.5)
    ax.legend(fontsize=7, loc="upper right", ncol=2)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, max(gdf["sig_rate"]) + 20)

plt.tight_layout()
fig.savefig(OUT_DIR / "domain_regression_summary.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: domain_regression_summary.png")

# ── Figure 2: Heatmap — Domain × Predictor (mean r, aggregated) ────────────

# Aggregate across all predictors: domain-level mean r and sig rate
all_labels = list(all_results.keys())

fig2, axes2 = plt.subplots(1, 2, figsize=(22, 8))
fig2.suptitle("Domain-level Summary: Mean r and Significance Rate\n"
              "(aggregated across all site-years within each domain)",
              fontsize=13, fontweight="bold")

for ax_i, (title, pred_subset) in enumerate([
    ("→ Richness", [l for l in all_labels if "Richness" in l]),
    ("→ LCBD", [l for l in all_labels if "LCBD" in l]),
]):
    ax = axes2[ax_i]
    domains_in = [d for d in domains_ordered
                  if any(d in all_results[p]["domain"].values for p in pred_subset if p in all_results)]

    r_mat = np.full((len(pred_subset), len(domains_in)), np.nan)
    sig_mat = np.full((len(pred_subset), len(domains_in)), np.nan)

    for i, pred in enumerate(pred_subset):
        if pred not in all_results:
            continue
        res = all_results[pred]
        for j, d in enumerate(domains_in):
            d_res = res[res["domain"] == d]
            if len(d_res) > 0:
                r_mat[i, j] = d_res["r"].mean()
                sig_mat[i, j] = (d_res["p_value"] < 0.05).mean() * 100

    im = ax.imshow(r_mat, aspect="auto", cmap="RdBu_r", vmin=-0.5, vmax=0.5)

    for i in range(len(pred_subset)):
        for j in range(len(domains_in)):
            if np.isnan(r_mat[i, j]):
                continue
            rv = r_mat[i, j]
            sr = sig_mat[i, j]
            color = "white" if abs(rv) > 0.3 else "black"
            ax.text(j, i, f"{rv:+.2f}\n({sr:.0f}%)", ha="center", va="center",
                    fontsize=7, color=color,
                    fontweight="bold" if sr > 30 else "normal")

    domain_labels = [f"{d}\n{DOMAIN_NAMES.get(d, '')}" for d in domains_in]
    ax.set_xticks(range(len(domains_in)))
    ax.set_xticklabels(domain_labels, fontsize=8)
    pred_labels = [p.split("→")[0].strip() for p in pred_subset]
    ax.set_yticks(range(len(pred_subset)))
    ax.set_yticklabels(pred_labels, fontsize=8)
    ax.set_title(title, fontsize=11, fontweight="bold")

    plt.colorbar(im, ax=ax, label="Mean Pearson r", shrink=0.8)

plt.tight_layout()
fig2.savefig(OUT_DIR / "domain_regression_heatmap.png", dpi=150, bbox_inches="tight")
print(f"Saved: domain_regression_heatmap.png")

# ── Print summary ───────────────────────────────────────────────────────────

print("\n" + "="*80)
print("DOMAIN-LEVEL SUMMARY")
print("="*80)

for label, res in sorted(all_results.items()):
    print(f"\n{label}:")
    for d in domains_ordered:
        d_res = res[res["domain"] == d]
        if len(d_res) == 0:
            continue
        sig = d_res[d_res["p_value"] < 0.05]
        sites = ", ".join(sorted(d_res["siteID"].unique()))
        print(f"  {d} ({DOMAIN_NAMES.get(d,''):15s}): {len(d_res):3d} site-yrs, "
              f"sig={len(sig):2d}/{len(d_res)} ({len(sig)/len(d_res):4.0%}), "
              f"mean r={d_res['r'].mean():+.3f}  [{sites}]")

print("\nDone.")
