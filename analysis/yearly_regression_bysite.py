"""
Within-Site Yearly Regression: Pooled Taxonomic ~ RS Diversity
===============================================================
For each site × year, run OLS regression.
Shows whether RS-taxonomy relationship holds WITHIN sites.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ── Load data ───────────────────────────────────────────────────────────────

alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness", "shannon", "simpson"]]

func = pd.read_csv("E:/neon_lidar/functional_diversity/functional_diversity_all.csv")
func = func[["siteID", "plotID", "year", "func_FRic", "func_FDiv", "func_FEve", "func_RaoQ"]]

fsd_plot = pd.read_csv("E:/neon_lidar/model_results/plot_level_complete.csv")
fsd_cols = ["siteID", "plotID", "fsd_year",
            "rumple_mean", "top_rugosity_mean", "mean_max_canopy_ht_mean",
            "FHD_mean", "LAI_mean", "vert_sd_mean", "GC_mean",
            "deepgap_fraction_mean", "VCI_mean"]
fsd_plot = fsd_plot[[c for c in fsd_cols if c in fsd_plot.columns]]
fsd_plot = fsd_plot.rename(columns={"fsd_year": "year"})

df_crown = pd.read_csv("E:/neon_lidar/spectral_diversity/deepforest_crown_diversity.csv")
df_crown_cols = ["siteID", "plotID", "year", "df_crown_RaoQ", "df_crown_FRic"]
df_crown = df_crown[[c for c in df_crown_cols if c in df_crown.columns]]

# Merge with pooled alpha
func_alpha = func.merge(alpha, on=["siteID", "plotID"], how="inner")
fsd_alpha = fsd_plot.merge(alpha, on=["siteID", "plotID"], how="inner")
dfcrown_alpha = df_crown.merge(alpha, on=["siteID", "plotID"], how="inner")

print(f"Func × Alpha: {len(func_alpha)} rows, {func_alpha.siteID.nunique()} sites")
print(f"FSD × Alpha: {len(fsd_alpha)} rows, {fsd_alpha.siteID.nunique()} sites")
print(f"Crown × Alpha: {len(dfcrown_alpha)} rows, {dfcrown_alpha.siteID.nunique()} sites")

# ── Site-year regression ────────────────────────────────────────────────────

def site_year_regression(df, x_col, y_col, min_n=8):
    """Run regression for each site × year."""
    results = []
    for (site, year), sub in df.groupby(["siteID", "year"]):
        sub = sub.dropna(subset=[x_col, y_col])
        if len(sub) < min_n:
            continue
        x, y = sub[x_col].values, sub[y_col].values
        if np.std(x) < 1e-10 or np.std(y) < 1e-10:
            continue
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        results.append({
            "siteID": site, "year": int(year), "n": len(sub),
            "r": r_value, "r2": r_value**2, "slope": slope,
            "p_value": p_value,
        })
    return pd.DataFrame(results)

# Key predictors
predictors = [
    ("func_RaoQ", "shannon", "Func RaoQ → Shannon", func_alpha),
    ("FHD_mean", "shannon", "FHD → Shannon", fsd_alpha),
    ("rumple_mean", "shannon", "Rumple → Shannon", fsd_alpha),
    ("vert_sd_mean", "shannon", "Vert SD → Shannon", fsd_alpha),
    ("func_RaoQ", "richness", "Func RaoQ → Richness", func_alpha),
    ("df_crown_RaoQ", "shannon", "Crown RaoQ (RGB) → Shannon", dfcrown_alpha),
]

all_results = {}
for x_col, y_col, label, df in predictors:
    if x_col not in df.columns:
        continue
    res = site_year_regression(df, x_col, y_col)
    if len(res) == 0:
        continue
    all_results[label] = res
    sig = res[res["p_value"] < 0.05]
    print(f"\n{label}: {len(res)} site-years, {len(sig)} significant ({len(sig)/len(res):.0%})")
    print(f"  Mean R²={res['r2'].mean():.3f}, mean r={res['r'].mean():.3f}")

# ── Figure 1: Heatmap (site × year) for top predictors ──────────────────────

from site_config import SITES

top_labels = ["Func RaoQ → Shannon", "FHD → Shannon", "Rumple → Shannon",
              "Func RaoQ → Richness"]
top_labels = [l for l in top_labels if l in all_results]

fig, axes = plt.subplots(len(top_labels), 1, figsize=(16, 4 * len(top_labels)))
if len(top_labels) == 1:
    axes = [axes]

fig.suptitle("Within-Site Yearly R²: Pooled Taxonomic ~ RS Diversity\n"
             "(each cell = one site-year regression; black border = p < 0.05)",
             fontsize=13, fontweight="bold", y=1.02)

for idx, label in enumerate(top_labels):
    ax = axes[idx]
    res = all_results[label]

    sites_in_data = sorted(res["siteID"].unique())
    years_in_data = sorted(res["year"].unique())

    # Build matrix
    r_matrix = np.full((len(sites_in_data), len(years_in_data)), np.nan)
    p_matrix = np.full((len(sites_in_data), len(years_in_data)), np.nan)
    n_matrix = np.full((len(sites_in_data), len(years_in_data)), np.nan)

    site_idx = {s: i for i, s in enumerate(sites_in_data)}
    year_idx = {y: i for i, y in enumerate(years_in_data)}

    for _, row in res.iterrows():
        si, yi = site_idx[row["siteID"]], year_idx[row["year"]]
        r_matrix[si, yi] = row["r"]
        p_matrix[si, yi] = row["p_value"]
        n_matrix[si, yi] = row["n"]

    # Plot heatmap using r (not R²) to show direction
    im = ax.imshow(r_matrix, aspect="auto", cmap="RdBu_r", vmin=-0.8, vmax=0.8)

    # Add text and borders for significant
    for i in range(len(sites_in_data)):
        for j in range(len(years_in_data)):
            if np.isnan(r_matrix[i, j]):
                continue
            r_val = r_matrix[i, j]
            n_val = int(n_matrix[i, j]) if not np.isnan(n_matrix[i, j]) else 0
            color = "white" if abs(r_val) > 0.4 else "black"
            ax.text(j, i, f"{r_val:.2f}\n({n_val})", ha="center", va="center",
                    fontsize=6, color=color)
            if p_matrix[i, j] < 0.05:
                rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, linewidth=2,
                                     edgecolor="black", facecolor="none")
                ax.add_patch(rect)

    ax.set_xticks(range(len(years_in_data)))
    ax.set_xticklabels(years_in_data, fontsize=8)
    ax.set_yticks(range(len(sites_in_data)))
    ax.set_yticklabels(sites_in_data, fontsize=9)
    ax.set_title(label, fontsize=11, fontweight="bold")

    plt.colorbar(im, ax=ax, label="Pearson r", shrink=0.8)

plt.tight_layout()
fig.savefig("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs/yearly_regression_bysite_heatmap.png",
            dpi=150, bbox_inches="tight")
print(f"\nHeatmap saved: docs/yearly_regression_bysite_heatmap.png")

# ── Figure 2: Site-level scatter facets for best predictor (Func RaoQ) ──────

if "Func RaoQ → Shannon" in all_results:
    res = all_results["Func RaoQ → Shannon"]
    sites_with_data = sorted(res["siteID"].unique())

    n_sites = len(sites_with_data)
    n_cols = 4
    n_rows = (n_sites + n_cols - 1) // n_cols

    fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
    fig2.suptitle("Within-Site: Func RaoQ → Shannon (each color = different year)\n"
                  "Taxonomic diversity pooled; structural diversity per year",
                  fontsize=14, fontweight="bold")
    axes2 = axes2.flatten()

    cmap = plt.cm.viridis
    all_years = sorted(func_alpha["year"].unique())
    year_colors = {y: cmap(i / max(len(all_years)-1, 1)) for i, y in enumerate(all_years)}

    for i, site in enumerate(sites_with_data):
        ax = axes2[i]
        site_data = func_alpha[func_alpha["siteID"] == site]

        for year in sorted(site_data["year"].unique()):
            yd = site_data[site_data["year"] == year].dropna(subset=["func_RaoQ", "shannon"])
            if len(yd) < 3:
                continue
            ax.scatter(yd["func_RaoQ"], yd["shannon"], c=[year_colors[year]],
                       s=25, alpha=0.7, label=str(int(year)))

            # Regression line if enough points
            if len(yd) >= 8:
                x, y = yd["func_RaoQ"].values, yd["shannon"].values
                slope, intercept, r, p, _ = linregress(x, y)
                xline = np.linspace(x.min(), x.max(), 50)
                ls = "-" if p < 0.05 else "--"
                ax.plot(xline, slope * xline + intercept, color=year_colors[year],
                        linewidth=1.5, linestyle=ls, alpha=0.7)

        ax.set_title(f"{site} ({len(site_data)} obs)", fontsize=10)
        ax.set_xlabel("Func RaoQ")
        ax.set_ylabel("Shannon")
        ax.legend(fontsize=6, ncol=2, loc="best")
        ax.grid(alpha=0.3)

    # Hide empty axes
    for j in range(i + 1, len(axes2)):
        axes2[j].set_visible(False)

    plt.tight_layout()
    fig2.savefig("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs/yearly_regression_bysite_scatter.png",
                 dpi=150, bbox_inches="tight")
    print(f"Scatter saved: docs/yearly_regression_bysite_scatter.png")

# ── Summary table ────────────────────────────────────────────────────────────

print("\n" + "="*70)
print("WITHIN-SITE SUMMARY (across all site-years)")
print("="*70)

all_res_combined = []
for label, res in all_results.items():
    sig = res[res["p_value"] < 0.05]
    pos_sig = sig[sig["r"] > 0]
    neg_sig = sig[sig["r"] < 0]
    print(f"\n{label}:")
    print(f"  Site-years: {len(res)}, Sig: {len(sig)} ({len(sig)/len(res):.0%})")
    print(f"  Positive sig: {len(pos_sig)}, Negative sig: {len(neg_sig)}")
    print(f"  Mean r: {res['r'].mean():.3f} ± {res['r'].std():.3f}")
    print(f"  Mean R²: {res['r2'].mean():.3f}")

    # Per-site summary
    for site in sorted(res["siteID"].unique()):
        sr = res[res["siteID"] == site]
        s_sig = sr[sr["p_value"] < 0.05]
        print(f"    {site}: {len(sr)} yrs, mean r={sr['r'].mean():.3f}, "
              f"sig={len(s_sig)}/{len(sr)}, mean n={sr['n'].mean():.0f}")

    res_copy = res.copy()
    res_copy["label"] = label
    all_res_combined.append(res_copy)

combined = pd.concat(all_res_combined, ignore_index=True)
combined.to_csv("E:/neon_lidar/model_results/yearly_regression_bysite.csv", index=False)
print(f"\nResults saved: {len(combined)} rows")
print("\nDone.")
