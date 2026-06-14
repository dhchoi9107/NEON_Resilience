"""
Domain & Site Heatmaps: Structural and Spectral separately
============================================================
8 heatmaps total:
  1) Richness ~ Structural (by domain)
  2) Richness ~ Structural (by site)
  3) Richness ~ Spectral (by domain)
  4) Richness ~ Spectral (by site)
  5) LCBD ~ Structural (by domain)
  6) LCBD ~ Structural (by site)
  7) LCBD ~ Spectral (by domain)
  8) LCBD ~ Spectral (by site)
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
    "D01": "Northeast", "D02": "Mid-Atlantic", "D03": "Southeast",
    "D05": "Great Lakes", "D07": "Appalachians", "D08": "Ozarks",
    "D10": "Rockies", "D16": "Pacific NW", "D17": "Pacific SW",
}
DOMAINS_ORDERED = ["D01", "D02", "D03", "D05", "D07", "D08", "D10", "D16", "D17"]

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

print(f"Structural: {len(fsd_taxo)} rows, Spectral: {len(vi_taxo)} rows")

# ── Predictor lists ─────────────────────────────────────────────────────────

struct_predictors = [
    ("FHD_mean", "FHD"),
    ("rumple_mean", "Rumple"),
    ("mean_max_canopy_ht_mean", "Canopy Ht"),
    ("max_canopy_ht_mean", "Max Ht"),
    ("vert_sd_mean", "Vert SD"),
    ("vertCV_mean", "Vert CV"),
    ("LAI_mean", "LAI"),
    ("GC_mean", "Gap Frac"),
    ("deepgap_fraction_mean", "Deep Gap"),
    ("VCI_mean", "VCI"),
    ("top_rugosity_mean", "Rugosity"),
    ("q95_mean", "Q95"),
    ("HeightRatio_mean", "Ht Ratio"),
]

spec_predictors = [
    ("NDVI_mean", "NDVI"),
    ("EVI_mean", "EVI"),
    ("ARVI_mean", "ARVI"),
    ("PRI_mean", "PRI"),
    ("SAVI_mean", "SAVI"),
    ("LAI_mean", "LAI(opt)"),
    ("fPAR_mean", "fPAR"),
    ("NDVI_sd", "NDVI SD"),
    ("EVI_sd", "EVI SD"),
    ("PRI_sd", "PRI SD"),
]

# ── Regression engine ───────────────────────────────────────────────────────

def run_regressions(df, x_col, y_col, group_col, year_col, min_n=8):
    """Run regression per group (aggregating across years within group)."""
    results = []
    for (grp, yr), sub in df.groupby([group_col, year_col]):
        sub = sub.dropna(subset=[x_col, y_col])
        if len(sub) < min_n:
            continue
        x, y = sub[x_col].values, sub[y_col].values
        if np.std(x) < 1e-10 or np.std(y) < 1e-10:
            continue
        _, _, r_value, p_value, _ = linregress(x, y)
        results.append({"group": grp, "year": int(yr), "n": len(sub),
                         "r": r_value, "p_value": p_value})
    return pd.DataFrame(results)


def aggregate_by_group(reg_df):
    """Aggregate site-year results into group-level summary."""
    agg = []
    for grp, gdf in reg_df.groupby("group"):
        sig = gdf[gdf["p_value"] < 0.05]
        pos_sig = sig[sig["r"] > 0]
        neg_sig = sig[sig["r"] < 0]
        agg.append({
            "group": grp,
            "n_site_years": len(gdf),
            "mean_r": gdf["r"].mean(),
            "sig_rate": len(sig) / len(gdf) if len(gdf) > 0 else 0,
            "pos_sig": len(pos_sig),
            "neg_sig": len(neg_sig),
        })
    return pd.DataFrame(agg)

# ── Heatmap drawing ─────────────────────────────────────────────────────────

def draw_domain_site_heatmap(results_dict, row_labels, title, filename,
                              vmin=-0.6, vmax=0.6):
    """
    results_dict: {predictor_label: DataFrame with 'group', 'mean_r', 'sig_rate'}
    row_labels: ordered list of groups (domains or sites)
    """
    pred_labels = list(results_dict.keys())
    n_pred = len(pred_labels)
    n_rows = len(row_labels)

    r_mat = np.full((n_rows, n_pred), np.nan)
    sig_mat = np.full((n_rows, n_pred), np.nan)

    row_idx = {r: i for i, r in enumerate(row_labels)}

    for j, pred in enumerate(pred_labels):
        agg = results_dict[pred]
        for _, row in agg.iterrows():
            if row["group"] in row_idx:
                i = row_idx[row["group"]]
                r_mat[i, j] = row["mean_r"]
                sig_mat[i, j] = row["sig_rate"] * 100

    fig_h = max(4, 0.55 * n_rows + 2)
    fig_w = max(8, 0.9 * n_pred + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(r_mat, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)

    for i in range(n_rows):
        for j in range(n_pred):
            if np.isnan(r_mat[i, j]):
                continue
            rv = r_mat[i, j]
            sr = sig_mat[i, j]
            color = "white" if abs(rv) > 0.35 else "black"
            fontw = "bold" if sr >= 30 else "normal"
            ax.text(j, i, f"{rv:+.2f}\n({sr:.0f}%)", ha="center", va="center",
                    fontsize=7, color=color, fontweight=fontw)
            if sr >= 50:
                rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, linewidth=2.5,
                                     edgecolor="black", facecolor="none")
                ax.add_patch(rect)

    ax.set_xticks(range(n_pred))
    ax.set_xticklabels(pred_labels, fontsize=9, rotation=45, ha="right")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

    plt.colorbar(im, ax=ax, label="Mean Pearson r", shrink=0.8, pad=0.02)
    plt.tight_layout()
    fig.savefig(OUT_DIR / filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {filename}")

# ── Compute all regressions ─────────────────────────────────────────────────

def compute_all(df, predictors, y_col, group_col, year_col):
    """Returns {pred_label: aggregated_df}"""
    out = {}
    for col, label in predictors:
        if col not in df.columns:
            continue
        reg = run_regressions(df, col, y_col, group_col, year_col)
        if len(reg) == 0:
            continue
        agg = aggregate_by_group(reg)
        out[label] = agg
    return out

# Year column
fsd_taxo_yr = fsd_taxo.rename(columns={"fsd_year": "year_col"})
vi_taxo_yr = vi_taxo.copy()
vi_taxo_yr["year_col"] = vi_taxo_yr["year"]

# ── 1) Richness ~ Structural ───────────────────────────────────────────────

# By domain
res = compute_all(fsd_taxo_yr, struct_predictors, "richness", "domain", "year_col")
domain_labels = [f"{d} {DOMAIN_NAMES.get(d,'')}" for d in DOMAINS_ORDERED
                 if any(d in r["group"].values for r in res.values())]
domain_keys = [d for d in DOMAINS_ORDERED
               if any(d in r["group"].values for r in res.values())]
# Remap labels
res_mapped = {}
for k, v in res.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda d: f"{d} {DOMAIN_NAMES.get(d,'')}")
    res_mapped[k] = v2
draw_domain_site_heatmap(res_mapped, domain_labels,
    "Richness ~ Structural (by Domain)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_rich_struct_domain.png")

# By site (ordered by domain)
sites_ordered = []
for d in DOMAINS_ORDERED:
    sites_in_d = sorted([s for s, dom in SITE_DOMAIN.items()
                         if dom == d and s in fsd_taxo["siteID"].unique()])
    sites_ordered.extend(sites_in_d)
site_labels = [f"{s} ({SITE_DOMAIN.get(s,'')})" for s in sites_ordered]

res_site = compute_all(fsd_taxo_yr, struct_predictors, "richness", "siteID", "year_col")
res_site_mapped = {}
for k, v in res_site.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda s: f"{s} ({SITE_DOMAIN.get(s,'')})")
    res_site_mapped[k] = v2
draw_domain_site_heatmap(res_site_mapped, site_labels,
    "Richness ~ Structural (by Site)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_rich_struct_site.png")

# ── 2) Richness ~ Spectral ─────────────────────────────────────────────────

res = compute_all(vi_taxo_yr, spec_predictors, "richness", "domain", "year_col")
res_mapped = {}
for k, v in res.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda d: f"{d} {DOMAIN_NAMES.get(d,'')}")
    res_mapped[k] = v2
domain_labels_vi = [f"{d} {DOMAIN_NAMES.get(d,'')}" for d in DOMAINS_ORDERED
                    if any(d in r["group"].values for r in res.values())]
draw_domain_site_heatmap(res_mapped, domain_labels_vi,
    "Richness ~ Spectral VI (by Domain)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_rich_spec_domain.png")

sites_vi = sorted(vi_taxo["siteID"].unique())
sites_vi_ordered = [s for s in sites_ordered if s in sites_vi]
# Add any missing
for s in sites_vi:
    if s not in sites_vi_ordered:
        sites_vi_ordered.append(s)
site_labels_vi = [f"{s} ({SITE_DOMAIN.get(s,'')})" for s in sites_vi_ordered]

res_site = compute_all(vi_taxo_yr, spec_predictors, "richness", "siteID", "year_col")
res_site_mapped = {}
for k, v in res_site.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda s: f"{s} ({SITE_DOMAIN.get(s,'')})")
    res_site_mapped[k] = v2
draw_domain_site_heatmap(res_site_mapped, site_labels_vi,
    "Richness ~ Spectral VI (by Site)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_rich_spec_site.png")

# ── 3) LCBD ~ Structural ───────────────────────────────────────────────────

res = compute_all(fsd_taxo_yr, struct_predictors, "lcbd_bray", "domain", "year_col")
res_mapped = {}
for k, v in res.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda d: f"{d} {DOMAIN_NAMES.get(d,'')}")
    res_mapped[k] = v2
draw_domain_site_heatmap(res_mapped, domain_labels,
    "LCBD ~ Structural (by Domain)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_lcbd_struct_domain.png")

res_site = compute_all(fsd_taxo_yr, struct_predictors, "lcbd_bray", "siteID", "year_col")
res_site_mapped = {}
for k, v in res_site.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda s: f"{s} ({SITE_DOMAIN.get(s,'')})")
    res_site_mapped[k] = v2
draw_domain_site_heatmap(res_site_mapped, site_labels,
    "LCBD ~ Structural (by Site)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_lcbd_struct_site.png")

# ── 4) LCBD ~ Spectral ─────────────────────────────────────────────────────

res = compute_all(vi_taxo_yr, spec_predictors, "lcbd_bray", "domain", "year_col")
res_mapped = {}
for k, v in res.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda d: f"{d} {DOMAIN_NAMES.get(d,'')}")
    res_mapped[k] = v2
draw_domain_site_heatmap(res_mapped, domain_labels_vi,
    "LCBD ~ Spectral VI (by Domain)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_lcbd_spec_domain.png")

res_site = compute_all(vi_taxo_yr, spec_predictors, "lcbd_bray", "siteID", "year_col")
res_site_mapped = {}
for k, v in res_site.items():
    v2 = v.copy()
    v2["group"] = v2["group"].map(lambda s: f"{s} ({SITE_DOMAIN.get(s,'')})")
    res_site_mapped[k] = v2
draw_domain_site_heatmap(res_site_mapped, site_labels_vi,
    "LCBD ~ Spectral VI (by Site)\ncell = mean r, (sig%); black border ≥ 50%",
    "heatmap_lcbd_spec_site.png")

print("\nAll 8 heatmaps saved.")
