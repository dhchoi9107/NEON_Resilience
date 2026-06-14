"""
Yearly Regression: Pooled Taxonomic Diversity ~ RS Diversity per Year
======================================================================
- Taxonomic diversity: pooled across all years (fixed per plot)
- RS-based diversity: structural (FSD) and spectral (crown) per year
- For each year, run OLS regression and report R², slope, significance

Output:
  - docs/yearly_regression.png (multi-panel figure)
  - E:/neon_lidar/model_results/yearly_regression_results.csv
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, linregress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ── Load data ───────────────────────────────────────────────────────────────

# 1) Pooled taxonomic alpha diversity (fixed per plot)
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness", "shannon", "simpson"]]
print(f"Pooled alpha: {len(alpha)} plots")

# 2) Structural diversity from FSD (per plot-year)
func = pd.read_csv("E:/neon_lidar/functional_diversity/functional_diversity_all.csv")
func = func[["siteID", "plotID", "year", "func_FRic", "func_FDiv", "func_FEve", "func_RaoQ"]]
print(f"Functional diversity: {len(func)} plot-years")

# 3) FSD structural metrics (per plot-year) — read from 1m plot-level rasters
#    Use plot_level_complete which has mean/sd of 21 FSD bands
fsd_plot = pd.read_csv("E:/neon_lidar/model_results/plot_level_complete.csv")
fsd_cols = ["siteID", "plotID", "fsd_year",
            "rumple_mean", "top_rugosity_mean", "mean_max_canopy_ht_mean",
            "FHD_mean", "LAI_mean", "vert_sd_mean", "GC_mean",
            "deepgap_fraction_mean", "VCI_mean"]
fsd_plot = fsd_plot[[c for c in fsd_cols if c in fsd_plot.columns]]
fsd_plot = fsd_plot.rename(columns={"fsd_year": "year"})
print(f"FSD plot-level: {len(fsd_plot)} plot-years")

# 4) DeepForest crown spectral diversity (per plot-year)
df_crown = pd.read_csv("E:/neon_lidar/spectral_diversity/deepforest_crown_diversity.csv")
df_crown_cols = ["siteID", "plotID", "year", "n_detections", "mean_score",
                 "df_crown_n", "df_crown_FRic", "df_crown_RaoQ", "df_crown_FDiv"]
df_crown = df_crown[[c for c in df_crown_cols if c in df_crown.columns]]
print(f"DeepForest crowns: {len(df_crown)} plot-years")

# 5) Hyperspectral crown spectral (fewer years but richer)
crown_spec = pd.read_csv("E:/neon_lidar/spectral_diversity/crown_spectral_diversity.csv")
crown_spec_cols = ["siteID", "plotID", "year", "crown_n", "crown_FRic", "crown_RaoQ", "crown_FDiv"]
crown_spec = crown_spec[[c for c in crown_spec_cols if c in crown_spec.columns]]
print(f"Hyperspectral crown: {len(crown_spec)} plot-years")

# ── Merge with pooled taxonomic ────────────────────────────────────────────

def merge_with_alpha(rs_df, year_col="year"):
    """Merge RS data with pooled alpha diversity."""
    merged = rs_df.merge(alpha, on=["siteID", "plotID"], how="inner")
    return merged

func_alpha = merge_with_alpha(func)
print(f"\nFunctional × Alpha: {len(func_alpha)} rows")

fsd_alpha = merge_with_alpha(fsd_plot)
print(f"FSD × Alpha: {len(fsd_alpha)} rows")

dfcrown_alpha = merge_with_alpha(df_crown)
print(f"DeepForest × Alpha: {len(dfcrown_alpha)} rows")

crown_alpha = merge_with_alpha(crown_spec)
print(f"HyperCrown × Alpha: {len(crown_alpha)} rows")

# ── Yearly regression function ─────────────────────────────────────────────

def yearly_regression(df, x_col, y_col, min_n=15):
    """Run regression for each year. Returns DataFrame of results."""
    results = []
    for year in sorted(df["year"].unique()):
        sub = df[df["year"] == year].dropna(subset=[x_col, y_col])
        if len(sub) < min_n:
            continue
        x, y = sub[x_col].values, sub[y_col].values
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        rho, p_spear = spearmanr(x, y)
        results.append({
            "year": int(year),
            "n": len(sub),
            "n_sites": sub["siteID"].nunique(),
            "r": r_value,
            "r2": r_value**2,
            "slope": slope,
            "p_value": p_value,
            "spearman_rho": rho,
            "spearman_p": p_spear,
            "x_col": x_col,
            "y_col": y_col,
        })
    return pd.DataFrame(results)

# ── Run all regressions ────────────────────────────────────────────────────

print("\n" + "="*70)
print("YEARLY REGRESSION: Pooled Taxonomic ~ RS Diversity")
print("="*70)

all_results = []

# A) Structural FSD metrics ~ Shannon
struct_pairs = [
    ("rumple_mean", "shannon", "Rumple", fsd_alpha),
    ("FHD_mean", "shannon", "FHD", fsd_alpha),
    ("LAI_mean", "shannon", "LAI", fsd_alpha),
    ("vert_sd_mean", "shannon", "Vert SD", fsd_alpha),
    ("mean_max_canopy_ht_mean", "shannon", "Canopy Ht", fsd_alpha),
    ("GC_mean", "shannon", "Gap fraction", fsd_alpha),
]

# B) Functional diversity ~ Shannon
func_pairs = [
    ("func_FRic", "shannon", "Func FRic", func_alpha),
    ("func_RaoQ", "shannon", "Func RaoQ", func_alpha),
    ("func_FDiv", "shannon", "Func FDiv", func_alpha),
]

# C) Crown spectral ~ Shannon
crown_pairs = [
    ("df_crown_RaoQ", "shannon", "Crown RaoQ (RGB)", dfcrown_alpha),
    ("df_crown_FRic", "shannon", "Crown FRic (RGB)", dfcrown_alpha),
]

# D) Hyperspectral crown ~ Shannon
hyper_pairs = [
    ("crown_RaoQ", "shannon", "Crown RaoQ (Hyper)", crown_alpha),
    ("crown_FRic", "shannon", "Crown FRic (Hyper)", crown_alpha),
]

# Also do richness as response
richness_pairs = [
    ("func_RaoQ", "richness", "Func RaoQ→Rich", func_alpha),
    ("FHD_mean", "richness", "FHD→Rich", fsd_alpha),
    ("rumple_mean", "richness", "Rumple→Rich", fsd_alpha),
]

all_pairs = struct_pairs + func_pairs + crown_pairs + hyper_pairs + richness_pairs

for x_col, y_col, label, df in all_pairs:
    if x_col not in df.columns or y_col not in df.columns:
        print(f"  SKIP {label}: column missing")
        continue
    res = yearly_regression(df, x_col, y_col)
    if len(res) == 0:
        print(f"  SKIP {label}: no valid years")
        continue
    res["label"] = label
    all_results.append(res)

    sig = res[res["p_value"] < 0.05]
    print(f"\n{label} ({x_col} ~ {y_col}):")
    print(f"  Years: {len(res)}, Sig: {len(sig)}/{len(res)}")
    print(f"  R² range: {res['r2'].min():.3f} - {res['r2'].max():.3f}, mean: {res['r2'].mean():.3f}")
    for _, row in res.iterrows():
        star = "*" if row["p_value"] < 0.05 else " "
        print(f"    {int(row['year'])}: n={int(row['n']):3d}, R²={row['r2']:.3f}, "
              f"r={row['r']:.3f}, ρ={row['spearman_rho']:.3f} {star}")

results_df = pd.concat(all_results, ignore_index=True)
results_df.to_csv("E:/neon_lidar/model_results/yearly_regression_results.csv", index=False)
print(f"\nResults saved: {len(results_df)} rows")

# ── Figure ──────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(3, 2, figsize=(14, 16))
fig.suptitle("Yearly Regression: Pooled Taxonomic Diversity ~ RS Diversity\n"
             "(taxonomic diversity pooled across all years; RS measured per year)",
             fontsize=13, fontweight="bold")

# Panel configs: (labels to plot, title, ax)
panel_configs = [
    # Top left: Structural → Shannon
    (["Rumple", "FHD", "Canopy Ht", "Vert SD"], "Structural Metrics → Shannon", axes[0, 0]),
    # Top right: Structural → Richness
    (["Rumple→Rich", "FHD→Rich"], "Structural Metrics → Richness", axes[0, 1]),
    # Mid left: Functional → Shannon
    (["Func FRic", "Func RaoQ", "Func FDiv"], "Functional Diversity → Shannon", axes[1, 0]),
    # Mid right: Func → Richness
    (["Func RaoQ→Rich"], "Functional Diversity → Richness", axes[1, 1]),
    # Bot left: Crown spectral → Shannon
    (["Crown RaoQ (RGB)", "Crown FRic (RGB)"], "Crown Spectral (RGB) → Shannon", axes[2, 0]),
    # Bot right: Hyperspectral crown → Shannon
    (["Crown RaoQ (Hyper)", "Crown FRic (Hyper)"], "Crown Spectral (Hyper) → Shannon", axes[2, 1]),
]

colors = plt.cm.tab10(np.linspace(0, 1, 10))

for labels_to_plot, title, ax in panel_configs:
    for i, label in enumerate(labels_to_plot):
        sub = results_df[results_df["label"] == label]
        if len(sub) == 0:
            continue
        c = colors[i % len(colors)]
        ax.plot(sub["year"], sub["r2"], "o-", color=c, label=label, markersize=5)
        # Mark significant years
        sig = sub[sub["p_value"] < 0.05]
        if len(sig) > 0:
            ax.scatter(sig["year"], sig["r2"], color=c, s=80, zorder=5,
                       edgecolors="black", linewidths=1.5)

    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Year")
    ax.set_ylabel("R²")
    ax.set_ylim(-0.02, max(0.3, ax.get_ylim()[1]))
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs/yearly_regression.png",
            dpi=150, bbox_inches="tight")
print(f"\nFigure saved: docs/yearly_regression.png")

# ── Summary: temporal stability ─────────────────────────────────────────────

print("\n" + "="*70)
print("TEMPORAL STABILITY SUMMARY")
print("="*70)

for label in results_df["label"].unique():
    sub = results_df[results_df["label"] == label]
    if len(sub) < 3:
        continue
    cv_r2 = sub["r2"].std() / (sub["r2"].mean() + 1e-6)
    sig_rate = (sub["p_value"] < 0.05).mean()
    print(f"  {label:25s}: mean R²={sub['r2'].mean():.3f} ± {sub['r2'].std():.3f}, "
          f"CV={cv_r2:.2f}, sig rate={sig_rate:.0%} ({len(sub)} years)")

# ── Scatter plot: best year per predictor ────────────────────────────────────

fig2, axes2 = plt.subplots(2, 3, figsize=(16, 10))
fig2.suptitle("Best-year Scatter: Pooled Taxonomic ~ RS Diversity", fontsize=13, fontweight="bold")

scatter_configs = [
    ("FHD_mean", "shannon", "FHD", fsd_alpha, axes2[0, 0]),
    ("rumple_mean", "shannon", "Rumple", fsd_alpha, axes2[0, 1]),
    ("func_RaoQ", "shannon", "Func RaoQ", func_alpha, axes2[0, 2]),
    ("func_FRic", "shannon", "Func FRic", func_alpha, axes2[1, 0]),
    ("df_crown_RaoQ", "shannon", "Crown RaoQ (RGB)", dfcrown_alpha, axes2[1, 1]),
    ("vert_sd_mean", "shannon", "Vert SD", fsd_alpha, axes2[1, 2]),
]

for x_col, y_col, label, df, ax in scatter_configs:
    if x_col not in df.columns:
        ax.set_visible(False)
        continue

    sub = results_df[results_df["label"] == label]
    if len(sub) == 0:
        ax.set_visible(False)
        continue
    best_year = sub.loc[sub["r2"].idxmax(), "year"]

    data = df[df["year"] == best_year].dropna(subset=[x_col, y_col])
    if len(data) < 10:
        ax.set_visible(False)
        continue

    # Color by site
    sites = sorted(data["siteID"].unique())
    site_colors = {s: plt.cm.tab20(i / max(len(sites)-1, 1)) for i, s in enumerate(sites)}

    for site in sites:
        s = data[data["siteID"] == site]
        ax.scatter(s[x_col], s[y_col], c=[site_colors[site]], s=20, alpha=0.6, label=site)

    # Regression line
    x, y = data[x_col].values, data[y_col].values
    slope, intercept, r, p, _ = linregress(x, y)
    xline = np.linspace(x.min(), x.max(), 100)
    ax.plot(xline, slope * xline + intercept, "r-", linewidth=2)

    ax.set_title(f"{label} ({int(best_year)})\nR²={r**2:.3f}, p={p:.1e}, n={len(data)}", fontsize=10)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.grid(alpha=0.3)

plt.tight_layout()
fig2.savefig("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs/yearly_regression_scatter.png",
             dpi=150, bbox_inches="tight")
print(f"Scatter figure saved: docs/yearly_regression_scatter.png")

print("\nDone.")
