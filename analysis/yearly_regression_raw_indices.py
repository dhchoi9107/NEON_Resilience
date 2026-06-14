"""
Within-Site Yearly Regression: Raw Indices
===========================================
Individual structural metrics (FHD, Height, Rumple, etc.) and
vegetation indices (NDVI, EVI, PRI, etc.) vs taxonomic diversity.

4 categories:
  1) Richness ~ Structural (individual LiDAR metrics)
  2) Richness ~ Spectral (individual VIs)
  3) LCBD ~ Structural
  4) LCBD ~ Spectral
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

OUT_DIR = Path("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs")

# ── Load data ───────────────────────────────────────────────────────────────

# Pooled taxonomic alpha (fixed per plot)
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness"]]

# LCBD (fixed per plot)
lcbd_src = pd.read_csv("E:/neon_lidar/model_results/plot_level_dbh10.csv",
                        usecols=["siteID", "plotID", "lcbd_bray"])
lcbd_src = lcbd_src.drop_duplicates("plotID")

# Merge alpha + LCBD
taxo = alpha.merge(lcbd_src, on=["siteID", "plotID"], how="inner")
print(f"Taxonomic (pooled): {len(taxo)} plots")

# ── Structural: FSD plot-level raw metrics (per plot-year) ──────────────────

fsd = pd.read_csv("E:/neon_lidar/model_results/plot_level_complete.csv")
struct_metrics = [
    ("FHD_mean",                "FHD"),
    ("rumple_mean",             "Rumple"),
    ("mean_max_canopy_ht_mean", "Mean Canopy Ht"),
    ("max_canopy_ht_mean",      "Max Canopy Ht"),
    ("vert_sd_mean",            "Vertical SD"),
    ("vertCV_mean",             "Vertical CV"),
    ("LAI_mean",                "LAI (LiDAR)"),
    ("GC_mean",                 "Gap Fraction"),
    ("deepgap_fraction_mean",   "Deep Gap Frac"),
    ("VCI_mean",                "VCI"),
    ("top_rugosity_mean",       "Top Rugosity"),
    ("meanH_mean",              "Mean Height"),
    ("q50_mean",                "Height Q50"),
    ("q95_mean",                "Height Q95"),
    ("HeightRatio_mean",        "Height Ratio"),
]

fsd_slim = fsd[["siteID", "plotID", "fsd_year"] +
               [m[0] for m in struct_metrics if m[0] in fsd.columns]].copy()
fsd_slim = fsd_slim.rename(columns={"fsd_year": "year"})
fsd_taxo = fsd_slim.merge(taxo, on=["siteID", "plotID"], how="inner")
print(f"Structural × Taxo: {len(fsd_taxo)} plot-years")

# ── Spectral: VI plot-level (per plot-year, 1m grain) ───────────────────────

vi = pd.read_csv("E:/neon_lidar/spectral_diversity/plot_spectral_1m.csv")
vi = vi[vi["grain_m"] == 1].copy()

spec_metrics = [
    ("NDVI_mean",  "NDVI"),
    ("EVI_mean",   "EVI"),
    ("ARVI_mean",  "ARVI"),
    ("PRI_mean",   "PRI"),
    ("SAVI_mean",  "SAVI"),
    ("LAI_mean",   "LAI (Optical)"),
    ("fPAR_mean",  "fPAR"),
    ("NDVI_sd",    "NDVI SD"),
    ("EVI_sd",     "EVI SD"),
    ("PRI_sd",     "PRI SD"),
]

vi_slim = vi[["siteID", "plotID", "year"] +
             [m[0] for m in spec_metrics if m[0] in vi.columns]].copy()
vi_taxo = vi_slim.merge(taxo, on=["siteID", "plotID"], how="inner")
print(f"Spectral × Taxo: {len(vi_taxo)} plot-years, years={sorted(vi_taxo.year.unique())}")

# ── Regression engine ───────────────────────────────────────────────────────

def site_year_regression(df, x_col, y_col, min_n=8):
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
            "r": r_value, "r2": r_value**2, "p_value": p_value,
        })
    return pd.DataFrame(results)

# ── Heatmap drawing ─────────────────────────────────────────────────────────

def draw_heatmap(res_df, title, ax):
    sites = sorted(res_df["siteID"].unique())
    years = sorted(res_df["year"].unique())
    r_mat = np.full((len(sites), len(years)), np.nan)
    p_mat = np.full((len(sites), len(years)), np.nan)
    si_map = {s: i for i, s in enumerate(sites)}
    yi_map = {y: i for i, y in enumerate(years)}
    for _, row in res_df.iterrows():
        si, yi = si_map[row["siteID"]], yi_map[row["year"]]
        r_mat[si, yi] = row["r"]
        p_mat[si, yi] = row["p_value"]
    im = ax.imshow(r_mat, aspect="auto", cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    for i in range(len(sites)):
        for j in range(len(years)):
            if np.isnan(r_mat[i, j]):
                continue
            rv = r_mat[i, j]
            color = "white" if abs(rv) > 0.45 else "black"
            ax.text(j, i, f"{rv:.2f}", ha="center", va="center",
                    fontsize=5.5, color=color,
                    fontweight="bold" if abs(rv) > 0.3 else "normal")
            if p_mat[i, j] < 0.05:
                rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, linewidth=2,
                                     edgecolor="black", facecolor="none")
                ax.add_patch(rect)
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, fontsize=7, rotation=45)
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites, fontsize=8)
    ax.set_title(title, fontsize=9, fontweight="bold")
    return im

# ── Run all and collect ─────────────────────────────────────────────────────

all_results = []

# Category 1: Richness ~ Structural
print("\n" + "="*70)
print("1) Richness ~ Structural (individual metrics)")
print("="*70)
cat1_results = {}
for col, label in struct_metrics:
    if col not in fsd_taxo.columns:
        continue
    res = site_year_regression(fsd_taxo, col, "richness")
    if len(res) == 0:
        continue
    cat1_results[label] = res
    sig = res[res["p_value"] < 0.05]
    pos = sig[sig["r"] > 0]
    neg = sig[sig["r"] < 0]
    print(f"  {label:20s}: {len(res):3d} site-yrs, sig={len(sig):2d}/{len(res)} ({len(sig)/len(res):4.0%}), "
          f"mean r={res['r'].mean():+.3f}, +{len(pos)}/-{len(neg)}")
    res_copy = res.copy()
    res_copy["label"] = label
    res_copy["category"] = "rich_struct"
    all_results.append(res_copy)

# Category 2: Richness ~ Spectral (individual VIs)
print("\n" + "="*70)
print("2) Richness ~ Spectral (individual VIs)")
print("="*70)
cat2_results = {}
for col, label in spec_metrics:
    if col not in vi_taxo.columns:
        continue
    res = site_year_regression(vi_taxo, col, "richness")
    if len(res) == 0:
        continue
    cat2_results[label] = res
    sig = res[res["p_value"] < 0.05]
    pos = sig[sig["r"] > 0]
    neg = sig[sig["r"] < 0]
    print(f"  {label:20s}: {len(res):3d} site-yrs, sig={len(sig):2d}/{len(res)} ({len(sig)/len(res):4.0%}), "
          f"mean r={res['r'].mean():+.3f}, +{len(pos)}/-{len(neg)}")
    res_copy = res.copy()
    res_copy["label"] = label
    res_copy["category"] = "rich_spec"
    all_results.append(res_copy)

# Category 3: LCBD ~ Structural
print("\n" + "="*70)
print("3) LCBD (Bray-Curtis) ~ Structural")
print("="*70)
cat3_results = {}
for col, label in struct_metrics:
    if col not in fsd_taxo.columns:
        continue
    res = site_year_regression(fsd_taxo, col, "lcbd_bray")
    if len(res) == 0:
        continue
    cat3_results[label] = res
    sig = res[res["p_value"] < 0.05]
    pos = sig[sig["r"] > 0]
    neg = sig[sig["r"] < 0]
    print(f"  {label:20s}: {len(res):3d} site-yrs, sig={len(sig):2d}/{len(res)} ({len(sig)/len(res):4.0%}), "
          f"mean r={res['r'].mean():+.3f}, +{len(pos)}/-{len(neg)}")
    res_copy = res.copy()
    res_copy["label"] = label
    res_copy["category"] = "lcbd_struct"
    all_results.append(res_copy)

# Category 4: LCBD ~ Spectral
print("\n" + "="*70)
print("4) LCBD (Bray-Curtis) ~ Spectral (individual VIs)")
print("="*70)
cat4_results = {}
for col, label in spec_metrics:
    if col not in vi_taxo.columns:
        continue
    res = site_year_regression(vi_taxo, col, "lcbd_bray")
    if len(res) == 0:
        continue
    cat4_results[label] = res
    sig = res[res["p_value"] < 0.05]
    pos = sig[sig["r"] > 0]
    neg = sig[sig["r"] < 0]
    print(f"  {label:20s}: {len(res):3d} site-yrs, sig={len(sig):2d}/{len(res)} ({len(sig)/len(res):4.0%}), "
          f"mean r={res['r'].mean():+.3f}, +{len(pos)}/-{len(neg)}")
    res_copy = res.copy()
    res_copy["label"] = label
    res_copy["category"] = "lcbd_spec"
    all_results.append(res_copy)

# Save all results
combined = pd.concat(all_results, ignore_index=True)
combined.to_csv("E:/neon_lidar/model_results/yearly_regression_raw_indices.csv", index=False)
print(f"\nResults saved: {len(combined)} rows")

# ── FIGURE 1: Richness ~ Structural heatmaps ──────────────────────────────

n = len(cat1_results)
fig, axes = plt.subplots(n, 1, figsize=(10, 3.5 * n))
if n == 1: axes = [axes]
fig.suptitle("Richness ~ Individual Structural Metrics (within-site, per year)\n"
             "black border = p<0.05", fontsize=13, fontweight="bold", y=1.01)
for i, (label, res) in enumerate(cat1_results.items()):
    im = draw_heatmap(res, label, axes[i])
cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar_ax, label="Pearson r")
plt.tight_layout(rect=[0, 0, 0.92, 0.98])
fig.savefig(OUT_DIR / "raw_richness_structural.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: raw_richness_structural.png ({n} panels)")
plt.close()

# ── FIGURE 2: Richness ~ Spectral heatmaps ────────────────────────────────

n = len(cat2_results)
fig, axes = plt.subplots(n, 1, figsize=(14, 3.5 * n))
if n == 1: axes = [axes]
fig.suptitle("Richness ~ Individual Vegetation Indices (within-site, per year)\n"
             "black border = p<0.05", fontsize=13, fontweight="bold", y=1.01)
for i, (label, res) in enumerate(cat2_results.items()):
    im = draw_heatmap(res, label, axes[i])
cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar_ax, label="Pearson r")
plt.tight_layout(rect=[0, 0, 0.92, 0.98])
fig.savefig(OUT_DIR / "raw_richness_spectral.png", dpi=150, bbox_inches="tight")
print(f"Saved: raw_richness_spectral.png ({n} panels)")
plt.close()

# ── FIGURE 3: LCBD ~ Structural heatmaps ──────────────────────────────────

n = len(cat3_results)
fig, axes = plt.subplots(n, 1, figsize=(10, 3.5 * n))
if n == 1: axes = [axes]
fig.suptitle("LCBD (Bray-Curtis) ~ Individual Structural Metrics (within-site, per year)\n"
             "black border = p<0.05", fontsize=13, fontweight="bold", y=1.01)
for i, (label, res) in enumerate(cat3_results.items()):
    im = draw_heatmap(res, label, axes[i])
cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar_ax, label="Pearson r")
plt.tight_layout(rect=[0, 0, 0.92, 0.98])
fig.savefig(OUT_DIR / "raw_lcbd_structural.png", dpi=150, bbox_inches="tight")
print(f"Saved: raw_lcbd_structural.png ({n} panels)")
plt.close()

# ── FIGURE 4: LCBD ~ Spectral heatmaps ────────────────────────────────────

n = len(cat4_results)
fig, axes = plt.subplots(n, 1, figsize=(14, 3.5 * n))
if n == 1: axes = [axes]
fig.suptitle("LCBD (Bray-Curtis) ~ Individual Vegetation Indices (within-site, per year)\n"
             "black border = p<0.05", fontsize=13, fontweight="bold", y=1.01)
for i, (label, res) in enumerate(cat4_results.items()):
    im = draw_heatmap(res, label, axes[i])
cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar_ax, label="Pearson r")
plt.tight_layout(rect=[0, 0, 0.92, 0.98])
fig.savefig(OUT_DIR / "raw_lcbd_spectral.png", dpi=150, bbox_inches="tight")
print(f"Saved: raw_lcbd_spectral.png ({n} panels)")
plt.close()

# ── FIGURE 5: Summary bar chart ───────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(18, 12))
fig.suptitle("Summary: Within-Site Significance Rate by Individual Index\n"
             "(% of site-years with p<0.05; numbers = mean r)",
             fontsize=14, fontweight="bold")

panels = [
    ("Richness ~ Structural", cat1_results, axes[0, 0], "#2196F3"),
    ("Richness ~ Spectral (VI)", cat2_results, axes[0, 1], "#FF9800"),
    ("LCBD ~ Structural", cat3_results, axes[1, 0], "#4CAF50"),
    ("LCBD ~ Spectral (VI)", cat4_results, axes[1, 1], "#9C27B0"),
]

for title, results_dict, ax, color in panels:
    if not results_dict:
        ax.set_visible(False)
        continue
    labels = list(results_dict.keys())
    sig_rates = []
    mean_rs = []
    for label in labels:
        res = results_dict[label]
        sig_rates.append((res["p_value"] < 0.05).mean() * 100)
        mean_rs.append(res["r"].mean())

    x = np.arange(len(labels))
    bars = ax.bar(x, sig_rates, color=color, edgecolor="black", linewidth=0.5, alpha=0.8)
    for i, (sr, mr) in enumerate(zip(sig_rates, mean_rs)):
        ax.text(i, sr + 0.8, f"r={mr:+.2f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Significant site-years (%)")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.axhline(y=5, color="red", linestyle="--", alpha=0.5, label="5% chance")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, max(sig_rates) + 10 if sig_rates else 30)

plt.tight_layout()
fig.savefig(OUT_DIR / "raw_indices_summary.png", dpi=150, bbox_inches="tight")
print(f"Saved: raw_indices_summary.png")

print("\nDone.")
