"""
Publication-quality Forest Plots & Interaction Plots
=====================================================
1. Forest plot: univariate mixed model coefficients (all predictors)
2. Interaction plot: domain-specific slopes for top predictors
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings("ignore")

OUT_DIR = Path("C:/Users/star1/Documents/GitHub/NEON_Resilience/docs")
RESULT_DIR = Path("E:/neon_lidar/model_results")

DOMAIN_NAMES = {
    "D01": "Northeast", "D02": "Mid-Atlantic", "D03": "Southeast",
    "D05": "Great Lakes", "D07": "Appalachians", "D08": "Ozarks",
    "D10": "Rockies", "D16": "Pacific NW", "D17": "Pacific SW",
}
DOMAINS_ORDERED = ["D01", "D02", "D03", "D05", "D07", "D08", "D10", "D16", "D17"]

# ── Load results ───────────────────────────────────────────────────────────

uni_df = pd.read_csv(RESULT_DIR / "mixed_model_univariate.csv")
int_df = pd.read_csv(RESULT_DIR / "mixed_model_interaction.csv")

print(f"Univariate: {len(uni_df)} models")
print(f"Interaction: {len(int_df)} models")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Forest plot — all predictors, Richness & LCBD side by side
# ══════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(14, 10), sharey=False)
fig.suptitle("Mixed Model: Standardized Coefficients\n"
             "response ~ z(predictor) + C(domain) + (1|site)",
             fontsize=14, fontweight="bold", y=0.98)

for ax_idx, resp in enumerate(["Richness", "LCBD"]):
    ax = axes[ax_idx]
    sub = uni_df[uni_df["resp_label"] == resp].copy()
    sub = sub.sort_values("beta")

    y_pos = np.arange(len(sub))

    # Color by category + significance
    colors = []
    for _, r in sub.iterrows():
        if r["p"] >= 0.05:
            colors.append("#cccccc")
        elif r["category"] == "Structural":
            colors.append("#2166ac")
        else:
            colors.append("#b2182b")

    ax.barh(y_pos, sub["beta"], xerr=1.96 * sub["se"],
            color=colors, edgecolor="black", linewidth=0.4,
            capsize=2, height=0.65, zorder=3)
    ax.axvline(x=0, color="black", linewidth=1, linestyle="-", zorder=2)

    # Y labels with significance stars
    ylabels = []
    for _, r in sub.iterrows():
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        cat = "S" if r["category"] == "Structural" else "V"
        ylabels.append(f"{r['pred_label']} ({cat}) {sig}")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Standardized coefficient (β ± 95% CI)", fontsize=10)
    ax.set_title(resp, fontsize=13, fontweight="bold")
    ax.grid(axis="x", alpha=0.2, zorder=1)

    # R²pred annotations on right side
    for i, (_, r) in enumerate(sub.iterrows()):
        if r["r2_predictor"] >= 0.01:
            ax.text(ax.get_xlim()[1] * 0.98, i,
                    f"R²={r['r2_predictor']:.3f}",
                    va="center", ha="right", fontsize=6.5,
                    style="italic", color="#333333")

legend_elements = [
    Patch(facecolor="#2166ac", edgecolor="black", label="Structural (p<0.05)"),
    Patch(facecolor="#b2182b", edgecolor="black", label="Spectral (p<0.05)"),
    Patch(facecolor="#cccccc", edgecolor="black", label="Not significant"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=3,
           fontsize=10, bbox_to_anchor=(0.5, -0.01),
           frameon=True, edgecolor="black")

plt.tight_layout(rect=[0, 0.04, 1, 0.95])
fig.savefig(OUT_DIR / "fig_forest_plot.png", dpi=200, bbox_inches="tight")
print("Saved: fig_forest_plot.png")
plt.close()


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Interaction heatmap — structural + spectral combined
# ══════════════════════════════════════════════════════════════════════════

for resp_label in ["Richness", "LCBD"]:
    int_sub = int_df[int_df["resp_label"] == resp_label].copy()
    if int_sub.empty:
        continue

    # Sort: structural first, then spectral, by abs mean slope
    int_sub["abs_mean_slope"] = int_sub[[f"slope_{d}" for d in DOMAINS_ORDERED]].abs().mean(axis=1)
    int_sub = int_sub.sort_values(["category", "abs_mean_slope"], ascending=[True, False])

    pred_labels = int_sub["pred_label"].tolist()
    n_pred = len(pred_labels)
    n_dom = len(DOMAINS_ORDERED)

    slope_mat = np.full((n_dom, n_pred), np.nan)
    for j, (_, r) in enumerate(int_sub.iterrows()):
        for i, d in enumerate(DOMAINS_ORDERED):
            col = f"slope_{d}"
            if col in r and not pd.isna(r[col]):
                slope_mat[i, j] = r[col]

    fig_h = max(5, 0.55 * n_dom + 2.5)
    fig_w = max(10, 0.75 * n_pred + 4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = np.nanmax(np.abs(slope_mat))
    im = ax.imshow(slope_mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    # Annotate cells
    for i in range(n_dom):
        for j in range(n_pred):
            if np.isnan(slope_mat[i, j]):
                ax.text(j, i, "—", ha="center", va="center", fontsize=7, color="#999999")
                continue
            sv = slope_mat[i, j]
            color = "white" if abs(sv) > vmax * 0.55 else "black"
            fontw = "bold" if abs(sv) > vmax * 0.35 else "normal"
            ax.text(j, i, f"{sv:+.2f}", ha="center", va="center",
                    fontsize=7.5, color=color, fontweight=fontw)

    # Domain labels
    domain_labels = [f"{d} {DOMAIN_NAMES.get(d, '')}" for d in DOMAINS_ORDERED]
    ax.set_yticks(range(n_dom))
    ax.set_yticklabels(domain_labels, fontsize=9)

    # Predictor labels with category and LRT significance
    xlabels = []
    for _, r in int_sub.iterrows():
        lrt_p = r.get("lrt_p", np.nan)
        if pd.isna(lrt_p):
            sig = ""
        elif lrt_p < 0.001:
            sig = "***"
        elif lrt_p < 0.01:
            sig = "**"
        elif lrt_p < 0.05:
            sig = "*"
        else:
            sig = ""
        cat = "S" if r["category"] == "Structural" else "V"
        xlabels.append(f"{r['pred_label']}({cat}){sig}")

    ax.set_xticks(range(n_pred))
    ax.set_xticklabels(xlabels, fontsize=8, rotation=45, ha="right")

    # Vertical line separating structural and spectral
    n_struct = sum(1 for _, r in int_sub.iterrows() if r["category"] == "Structural")
    if 0 < n_struct < n_pred:
        ax.axvline(x=n_struct - 0.5, color="black", linewidth=1.5, linestyle="--")

    ax.set_title(f"{resp_label} ~ Predictor × Domain\n"
                 f"Domain-specific slopes from interaction model; "
                 f"S=Structural, V=Spectral; *=LRT p<0.05",
                 fontsize=11, fontweight="bold", pad=12)

    cbar = plt.colorbar(im, ax=ax, label="Domain-specific slope (β)",
                        shrink=0.85, pad=0.02)
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    fname = f"fig_interaction_heatmap_{resp_label.lower()}.png"
    fig.savefig(OUT_DIR / fname, dpi=200, bbox_inches="tight")
    print(f"Saved: {fname}")
    plt.close()


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3: Interaction slope plots — regression lines per domain
# ══════════════════════════════════════════════════════════════════════════

from site_config import SITE_DOMAIN
from scipy.stats import linregress

# Load actual data for plotting
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness"]]
lcbd_src = pd.read_csv("E:/neon_lidar/model_results/plot_level_dbh10.csv",
                        usecols=["siteID", "plotID", "lcbd_bray"])
lcbd_src = lcbd_src.drop_duplicates("plotID")
taxo = alpha.merge(lcbd_src, on=["siteID", "plotID"], how="inner")
taxo["domain"] = taxo["siteID"].map(SITE_DOMAIN)

fsd = pd.read_csv("E:/neon_lidar/model_results/plot_level_complete.csv")
fsd["domain"] = fsd["siteID"].map(SITE_DOMAIN)
fsd_taxo = fsd.merge(taxo[["siteID", "plotID", "richness", "lcbd_bray"]],
                      on=["siteID", "plotID"], how="inner", suffixes=("_fsd", ""))

vi = pd.read_csv("E:/neon_lidar/spectral_diversity/plot_spectral_1m.csv")
vi = vi[vi["grain_m"] == 1].copy()
vi["domain"] = vi["siteID"].map(SITE_DOMAIN)
vi_taxo = vi.merge(taxo[["siteID", "plotID", "richness", "lcbd_bray"]],
                    on=["siteID", "plotID"], how="inner")

domain_colors = {
    "D01": "#e41a1c", "D02": "#377eb8", "D03": "#4daf4a",
    "D05": "#984ea3", "D07": "#ff7f00", "D08": "#a65628",
    "D10": "#f781bf", "D16": "#999999", "D17": "#66c2a5",
}

# Top predictors for interaction slope plots
slope_plots = [
    ("Richness", "Rumple", "rumple_mean", fsd_taxo, "richness"),
    ("Richness", "Deep Gap", "deepgap_fraction_mean", fsd_taxo, "richness"),
    ("Richness", "LAI", "LAI_mean", fsd_taxo, "richness"),
    ("Richness", "NDVI", "NDVI_mean", vi_taxo, "richness"),
    ("Richness", "EVI", "EVI_mean", vi_taxo, "richness"),
    ("Richness", "ARVI", "ARVI_mean", vi_taxo, "richness"),
    ("LCBD", "Vert CV", "vertCV_mean", fsd_taxo, "lcbd_bray"),
    ("LCBD", "Gap Frac", "GC_mean", fsd_taxo, "lcbd_bray"),
    ("LCBD", "ARVI", "ARVI_mean", vi_taxo, "lcbd_bray"),
    ("LCBD", "NDVI", "NDVI_mean", vi_taxo, "lcbd_bray"),
]

# Richness interaction plots
rich_plots = [p for p in slope_plots if p[0] == "Richness"]
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle("Richness ~ Predictor: Domain-specific Regression Lines\n"
             "Each line = one domain; points = plot-year observations",
             fontsize=13, fontweight="bold", y=0.99)

for idx, (resp_label, pred_name, col, df, ycol) in enumerate(rich_plots):
    ax = axes.flatten()[idx]
    sub = df[[col, ycol, "domain"]].dropna()

    for d in DOMAINS_ORDERED:
        dsub = sub[sub["domain"] == d]
        if len(dsub) < 10:
            continue

        x, y = dsub[col].values, dsub[ycol].values
        ax.scatter(x, y, c=domain_colors[d], alpha=0.15, s=8, zorder=2)

        # Regression line
        slope, intercept, r, p, _ = linregress(x, y)
        x_line = np.linspace(x.min(), x.max(), 50)
        lw = 2.5 if abs(r) > 0.3 else 1.2
        ls = "-" if p < 0.05 else "--"
        ax.plot(x_line, intercept + slope * x_line,
                c=domain_colors[d], linewidth=lw, linestyle=ls,
                label=f"{d} (r={r:+.2f})", zorder=3)

    ax.set_xlabel(pred_name, fontsize=10)
    ax.set_ylabel("Richness", fontsize=10)
    ax.set_title(f"Richness ~ {pred_name}", fontsize=11, fontweight="bold")
    ax.legend(fontsize=6, loc="best", ncol=2, framealpha=0.8)
    ax.grid(alpha=0.2)

# Hide unused subplot
for idx in range(len(rich_plots), 6):
    axes.flatten()[idx].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_DIR / "fig_interaction_slopes_richness.png", dpi=200, bbox_inches="tight")
print("Saved: fig_interaction_slopes_richness.png")
plt.close()

# LCBD interaction plots
lcbd_plots = [p for p in slope_plots if p[0] == "LCBD"]
fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle("LCBD ~ Predictor: Domain-specific Regression Lines\n"
             "Each line = one domain; points = plot-year observations",
             fontsize=13, fontweight="bold", y=0.99)

for idx, (resp_label, pred_name, col, df, ycol) in enumerate(lcbd_plots):
    ax = axes.flatten()[idx]
    sub = df[[col, ycol, "domain"]].dropna()

    for d in DOMAINS_ORDERED:
        dsub = sub[sub["domain"] == d]
        if len(dsub) < 10:
            continue

        x, y = dsub[col].values, dsub[ycol].values
        ax.scatter(x, y, c=domain_colors[d], alpha=0.15, s=8, zorder=2)

        slope, intercept, r, p, _ = linregress(x, y)
        x_line = np.linspace(x.min(), x.max(), 50)
        lw = 2.5 if abs(r) > 0.3 else 1.2
        ls = "-" if p < 0.05 else "--"
        ax.plot(x_line, intercept + slope * x_line,
                c=domain_colors[d], linewidth=lw, linestyle=ls,
                label=f"{d} (r={r:+.2f})", zorder=3)

    ax.set_xlabel(pred_name, fontsize=10)
    ax.set_ylabel("LCBD (Bray-Curtis)", fontsize=10)
    ax.set_title(f"LCBD ~ {pred_name}", fontsize=11, fontweight="bold")
    ax.legend(fontsize=6, loc="best", ncol=2, framealpha=0.8)
    ax.grid(alpha=0.2)

plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_DIR / "fig_interaction_slopes_lcbd.png", dpi=200, bbox_inches="tight")
print("Saved: fig_interaction_slopes_lcbd.png")
plt.close()


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4: Variance partitioning — clean version
# ══════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(16, 8))
fig.suptitle("Variance Partitioning\n"
             "response ~ z(predictor) + C(domain) + (1|site)",
             fontsize=13, fontweight="bold", y=0.99)

for ax_idx, resp in enumerate(["Richness", "LCBD"]):
    ax = axes[ax_idx]
    vp = uni_df[uni_df["resp_label"] == resp].copy()
    vp = vp.sort_values("r2_predictor", ascending=True)

    y = np.arange(len(vp))
    labels = [f"{r['pred_label']} ({r['category'][:3]})" for _, r in vp.iterrows()]

    # Clamp negative r2_domain to 0 for display
    r2_pred = vp["r2_predictor"].values.clip(0)
    r2_dom = vp["r2_domain"].values.clip(0)
    r2_site = (vp["r2_conditional"] - vp["r2_marginal"]).values.clip(0)
    r2_resid = (1.0 - r2_pred - r2_dom - r2_site).clip(0)

    ax.barh(y, r2_pred, color="#d73027", edgecolor="black", linewidth=0.3,
            height=0.7, label="Predictor", zorder=3)
    ax.barh(y, r2_dom, left=r2_pred, color="#fc8d59",
            edgecolor="black", linewidth=0.3, height=0.7,
            label="Domain (fixed)", zorder=3)
    ax.barh(y, r2_site, left=r2_pred + r2_dom, color="#4393c3",
            edgecolor="black", linewidth=0.3, height=0.7,
            label="Site (random)", zorder=3)
    ax.barh(y, r2_resid, left=r2_pred + r2_dom + r2_site,
            color="#e0e0e0", edgecolor="black", linewidth=0.3,
            height=0.7, label="Residual", zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Proportion of variance", fontsize=10)
    ax.set_title(resp, fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
    ax.grid(axis="x", alpha=0.2, zorder=1)
    ax.set_xlim(0, 1.0)

plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_DIR / "fig_variance_partition.png", dpi=200, bbox_inches="tight")
print("Saved: fig_variance_partition.png")
plt.close()

print("\nAll figures saved.")
