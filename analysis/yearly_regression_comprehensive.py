"""
Comprehensive Within-Site Yearly Regression
=============================================
4 categories × multiple RS predictors → heatmap per predictor

Categories:
  1) Alpha (Richness) ~ Structural
  2) Alpha (Richness) ~ Spectral
  3) Beta (LCBD Bray) ~ Structural
  4) Beta (LCBD Bray) ~ Spectral

Each heatmap: sites (rows) × years (columns), cell = Pearson r
Black border = p < 0.05
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

# Pooled taxonomic (fixed per plot)
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness"]]
print(f"Alpha pooled: {len(alpha)} plots")

# LCBD (fixed per plot, from dbh10)
lcbd = pd.read_csv("E:/neon_lidar/model_results/plot_level_dbh10.csv",
                    usecols=["siteID", "plotID", "lcbd_bray", "lcbd_turnover", "lcbd_nestedness"])
lcbd = lcbd.drop_duplicates("plotID")
print(f"LCBD pooled: {len(lcbd)} plots")

# Structural: functional diversity (FSD-based, per plot-year) — 14 years
func = pd.read_csv("E:/neon_lidar/functional_diversity/functional_diversity_all.csv")
func = func[["siteID", "plotID", "year", "func_FRic", "func_FDiv", "func_FEve", "func_RaoQ"]]
print(f"Functional (structural): {len(func)} plot-years, years={sorted(func.year.unique())}")

# Spectral: DeepForest crown (RGB, per plot-year) — 14 years
df_crown = pd.read_csv("E:/neon_lidar/spectral_diversity/deepforest_crown_diversity.csv")
df_crown = df_crown[["siteID", "plotID", "year", "df_crown_n", "df_crown_FRic",
                      "df_crown_RaoQ", "df_crown_FDiv"]].copy()
# Drop rows where crown metrics are 0 or NaN (no valid crowns)
df_crown = df_crown[df_crown["df_crown_n"] > 0]
print(f"Crown spectral (RGB): {len(df_crown)} plot-years")

# Hyperspectral crown (fewer years but richer spectral info)
hyper_crown = pd.read_csv("E:/neon_lidar/spectral_diversity/crown_spectral_diversity.csv")
hyper_crown = hyper_crown[["siteID", "plotID", "year", "crown_n", "crown_FRic",
                            "crown_RaoQ", "crown_FDiv"]].copy()
hyper_crown = hyper_crown[hyper_crown["crown_n"] > 0]
print(f"Crown spectral (Hyper): {len(hyper_crown)} plot-years")

# ── Merge datasets ──────────────────────────────────────────────────────────

# Structural × Alpha
struct_alpha = func.merge(alpha, on=["siteID", "plotID"], how="inner")
# Structural × Beta
struct_beta = func.merge(lcbd, on=["siteID", "plotID"], how="inner")
# Spectral (RGB) × Alpha
spec_alpha = df_crown.merge(alpha, on=["siteID", "plotID"], how="inner")
# Spectral (RGB) × Beta
spec_beta = df_crown.merge(lcbd, on=["siteID", "plotID"], how="inner")
# Hyper × Alpha
hyper_alpha = hyper_crown.merge(alpha, on=["siteID", "plotID"], how="inner")
# Hyper × Beta
hyper_beta = hyper_crown.merge(lcbd, on=["siteID", "plotID"], how="inner")

print(f"\nStruct × Alpha: {len(struct_alpha)}")
print(f"Struct × Beta:  {len(struct_beta)}")
print(f"Spec × Alpha:   {len(spec_alpha)}")
print(f"Spec × Beta:    {len(spec_beta)}")
print(f"Hyper × Alpha:  {len(hyper_alpha)}")
print(f"Hyper × Beta:   {len(hyper_beta)}")

# ── Regression function ────────────────────────────────────────────────────

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

# ── Define all predictor-response pairs ─────────────────────────────────────

pairs = [
    # Category 1: Alpha (Richness) ~ Structural
    ("func_RaoQ",  "richness",  "Struct RaoQ → Richness",       struct_alpha, "alpha_struct"),
    ("func_FRic",  "richness",  "Struct FRic → Richness",       struct_alpha, "alpha_struct"),
    ("func_FDiv",  "richness",  "Struct FDiv → Richness",       struct_alpha, "alpha_struct"),
    ("func_FEve",  "richness",  "Struct FEve → Richness",       struct_alpha, "alpha_struct"),

    # Category 2: Alpha (Richness) ~ Spectral (RGB)
    ("df_crown_RaoQ", "richness", "Spec RaoQ (RGB) → Richness", spec_alpha, "alpha_spec"),
    ("df_crown_FRic", "richness", "Spec FRic (RGB) → Richness", spec_alpha, "alpha_spec"),
    ("df_crown_FDiv", "richness", "Spec FDiv (RGB) → Richness", spec_alpha, "alpha_spec"),

    # Category 2b: Alpha (Richness) ~ Spectral (Hyperspectral)
    ("crown_RaoQ", "richness", "Spec RaoQ (Hyper) → Richness",  hyper_alpha, "alpha_hyper"),
    ("crown_FRic", "richness", "Spec FRic (Hyper) → Richness",  hyper_alpha, "alpha_hyper"),

    # Category 3: Beta (LCBD) ~ Structural
    ("func_RaoQ",  "lcbd_bray", "Struct RaoQ → LCBD",           struct_beta, "beta_struct"),
    ("func_FRic",  "lcbd_bray", "Struct FRic → LCBD",           struct_beta, "beta_struct"),
    ("func_FDiv",  "lcbd_bray", "Struct FDiv → LCBD",           struct_beta, "beta_struct"),
    ("func_FEve",  "lcbd_bray", "Struct FEve → LCBD",           struct_beta, "beta_struct"),

    # Category 4: Beta (LCBD) ~ Spectral (RGB)
    ("df_crown_RaoQ", "lcbd_bray", "Spec RaoQ (RGB) → LCBD",   spec_beta, "beta_spec"),
    ("df_crown_FRic", "lcbd_bray", "Spec FRic (RGB) → LCBD",   spec_beta, "beta_spec"),
    ("df_crown_FDiv", "lcbd_bray", "Spec FDiv (RGB) → LCBD",   spec_beta, "beta_spec"),

    # Category 4b: Beta (LCBD) ~ Spectral (Hyperspectral)
    ("crown_RaoQ", "lcbd_bray", "Spec RaoQ (Hyper) → LCBD",    hyper_beta, "beta_hyper"),
    ("crown_FRic", "lcbd_bray", "Spec FRic (Hyper) → LCBD",    hyper_beta, "beta_hyper"),
]

# ── Run all regressions ────────────────────────────────────────────────────

all_results = {}
for x_col, y_col, label, df, cat in pairs:
    if x_col not in df.columns or y_col not in df.columns:
        print(f"  SKIP {label}: column missing")
        continue
    res = site_year_regression(df, x_col, y_col)
    if len(res) == 0:
        print(f"  SKIP {label}: no valid site-years")
        continue
    all_results[label] = (res, cat)
    sig = res[res["p_value"] < 0.05]
    print(f"{label}: {len(res)} site-yrs, sig={len(sig)}/{len(res)} ({len(sig)/len(res):.0%}), "
          f"mean r={res['r'].mean():.3f}")

# ── Draw heatmaps ──────────────────────────────────────────────────────────

def draw_heatmap(res_df, title, ax, vmin=-0.8, vmax=0.8):
    """Draw site × year heatmap on given axes."""
    sites = sorted(res_df["siteID"].unique())
    years = sorted(res_df["year"].unique())

    r_matrix = np.full((len(sites), len(years)), np.nan)
    p_matrix = np.full((len(sites), len(years)), np.nan)
    n_matrix = np.full((len(sites), len(years)), np.nan)

    site_idx = {s: i for i, s in enumerate(sites)}
    year_idx = {y: i for i, y in enumerate(years)}

    for _, row in res_df.iterrows():
        si, yi = site_idx[row["siteID"]], year_idx[row["year"]]
        r_matrix[si, yi] = row["r"]
        p_matrix[si, yi] = row["p_value"]
        n_matrix[si, yi] = row["n"]

    im = ax.imshow(r_matrix, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)

    for i in range(len(sites)):
        for j in range(len(years)):
            if np.isnan(r_matrix[i, j]):
                continue
            r_val = r_matrix[i, j]
            n_val = int(n_matrix[i, j])
            color = "white" if abs(r_val) > 0.45 else "black"
            ax.text(j, i, f"{r_val:.2f}", ha="center", va="center",
                    fontsize=5.5, color=color, fontweight="bold" if abs(r_val) > 0.3 else "normal")
            if p_matrix[i, j] < 0.05:
                rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, linewidth=2,
                                     edgecolor="black", facecolor="none")
                ax.add_patch(rect)

    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, fontsize=7, rotation=45)
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites, fontsize=8)
    ax.set_title(title, fontsize=10, fontweight="bold")

    return im

# ── Figure 1: Alpha (Richness) ~ Structural (4 panels) ─────────────────────

alpha_struct_labels = [l for l, (_, c) in all_results.items() if c == "alpha_struct"]
n_panels = len(alpha_struct_labels)
if n_panels > 0:
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))
    if n_panels == 1: axes = [axes]
    fig.suptitle("Alpha (Richness) ~ Structural Diversity\n"
                 "(within-site, per year; black border = p<0.05)",
                 fontsize=14, fontweight="bold", y=1.01)
    for i, label in enumerate(alpha_struct_labels):
        res, _ = all_results[label]
        im = draw_heatmap(res, label, axes[i])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    plt.tight_layout(rect=[0, 0, 0.91, 0.97])
    fig.savefig(OUT_DIR / "comprehensive_alpha_structural.png", dpi=150, bbox_inches="tight")
    print(f"\nSaved: comprehensive_alpha_structural.png ({n_panels} panels)")

# ── Figure 2: Alpha (Richness) ~ Spectral RGB (3 panels) ───────────────────

alpha_spec_labels = [l for l, (_, c) in all_results.items() if c == "alpha_spec"]
n_panels = len(alpha_spec_labels)
if n_panels > 0:
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))
    if n_panels == 1: axes = [axes]
    fig.suptitle("Alpha (Richness) ~ Spectral Diversity (RGB Crown)\n"
                 "(within-site, per year; black border = p<0.05)",
                 fontsize=14, fontweight="bold", y=1.01)
    for i, label in enumerate(alpha_spec_labels):
        res, _ = all_results[label]
        im = draw_heatmap(res, label, axes[i])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    plt.tight_layout(rect=[0, 0, 0.91, 0.97])
    fig.savefig(OUT_DIR / "comprehensive_alpha_spectral_rgb.png", dpi=150, bbox_inches="tight")
    print(f"Saved: comprehensive_alpha_spectral_rgb.png ({n_panels} panels)")

# ── Figure 2b: Alpha ~ Spectral Hyper ──────────────────────────────────────

alpha_hyper_labels = [l for l, (_, c) in all_results.items() if c == "alpha_hyper"]
n_panels = len(alpha_hyper_labels)
if n_panels > 0:
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))
    if n_panels == 1: axes = [axes]
    fig.suptitle("Alpha (Richness) ~ Spectral Diversity (Hyperspectral Crown)\n"
                 "(within-site, per year; black border = p<0.05)",
                 fontsize=14, fontweight="bold", y=1.01)
    for i, label in enumerate(alpha_hyper_labels):
        res, _ = all_results[label]
        im = draw_heatmap(res, label, axes[i])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    plt.tight_layout(rect=[0, 0, 0.91, 0.97])
    fig.savefig(OUT_DIR / "comprehensive_alpha_spectral_hyper.png", dpi=150, bbox_inches="tight")
    print(f"Saved: comprehensive_alpha_spectral_hyper.png ({n_panels} panels)")

# ── Figure 3: Beta (LCBD) ~ Structural (4 panels) ─────────────────────────

beta_struct_labels = [l for l, (_, c) in all_results.items() if c == "beta_struct"]
n_panels = len(beta_struct_labels)
if n_panels > 0:
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))
    if n_panels == 1: axes = [axes]
    fig.suptitle("Beta (LCBD Bray-Curtis) ~ Structural Diversity\n"
                 "(within-site, per year; black border = p<0.05)",
                 fontsize=14, fontweight="bold", y=1.01)
    for i, label in enumerate(beta_struct_labels):
        res, _ = all_results[label]
        im = draw_heatmap(res, label, axes[i])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    plt.tight_layout(rect=[0, 0, 0.91, 0.97])
    fig.savefig(OUT_DIR / "comprehensive_beta_structural.png", dpi=150, bbox_inches="tight")
    print(f"Saved: comprehensive_beta_structural.png ({n_panels} panels)")

# ── Figure 4: Beta (LCBD) ~ Spectral RGB (3 panels) ───────────────────────

beta_spec_labels = [l for l, (_, c) in all_results.items() if c == "beta_spec"]
n_panels = len(beta_spec_labels)
if n_panels > 0:
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))
    if n_panels == 1: axes = [axes]
    fig.suptitle("Beta (LCBD Bray-Curtis) ~ Spectral Diversity (RGB Crown)\n"
                 "(within-site, per year; black border = p<0.05)",
                 fontsize=14, fontweight="bold", y=1.01)
    for i, label in enumerate(beta_spec_labels):
        res, _ = all_results[label]
        im = draw_heatmap(res, label, axes[i])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    plt.tight_layout(rect=[0, 0, 0.91, 0.97])
    fig.savefig(OUT_DIR / "comprehensive_beta_spectral_rgb.png", dpi=150, bbox_inches="tight")
    print(f"Saved: comprehensive_beta_spectral_rgb.png ({n_panels} panels)")

# ── Figure 4b: Beta ~ Spectral Hyper ──────────────────────────────────────

beta_hyper_labels = [l for l, (_, c) in all_results.items() if c == "beta_hyper"]
n_panels = len(beta_hyper_labels)
if n_panels > 0:
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4.5 * n_panels))
    if n_panels == 1: axes = [axes]
    fig.suptitle("Beta (LCBD Bray-Curtis) ~ Spectral Diversity (Hyperspectral Crown)\n"
                 "(within-site, per year; black border = p<0.05)",
                 fontsize=14, fontweight="bold", y=1.01)
    for i, label in enumerate(beta_hyper_labels):
        res, _ = all_results[label]
        im = draw_heatmap(res, label, axes[i])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    plt.tight_layout(rect=[0, 0, 0.91, 0.97])
    fig.savefig(OUT_DIR / "comprehensive_beta_spectral_hyper.png", dpi=150, bbox_inches="tight")
    print(f"Saved: comprehensive_beta_spectral_hyper.png ({n_panels} panels)")

# ── Figure 5: Summary comparison (all categories in one) ───────────────────

fig, ax = plt.subplots(figsize=(14, 8))

categories = {
    "Alpha~Struct": "alpha_struct",
    "Alpha~Spec(RGB)": "alpha_spec",
    "Alpha~Spec(Hyper)": "alpha_hyper",
    "Beta~Struct": "beta_struct",
    "Beta~Spec(RGB)": "beta_spec",
    "Beta~Spec(Hyper)": "beta_hyper",
}

summary_data = []
for cat_label, cat_key in categories.items():
    labels_in_cat = [l for l, (_, c) in all_results.items() if c == cat_key]
    for label in labels_in_cat:
        res, _ = all_results[label]
        short = label.split("→")[0].strip()
        sig_rate = (res["p_value"] < 0.05).mean()
        mean_r = res["r"].mean()
        mean_r2 = res["r2"].mean()
        n_sy = len(res)
        pos_sig = ((res["p_value"] < 0.05) & (res["r"] > 0)).sum()
        neg_sig = ((res["p_value"] < 0.05) & (res["r"] < 0)).sum()
        summary_data.append({
            "category": cat_label, "predictor": short, "label": label,
            "n_site_years": n_sy, "mean_r": mean_r, "mean_r2": mean_r2,
            "sig_rate": sig_rate, "pos_sig": pos_sig, "neg_sig": neg_sig,
        })

summary_df = pd.DataFrame(summary_data)

# Bar chart: sig_rate by predictor, grouped by category
x_labels = summary_df["label"].values
x = np.arange(len(x_labels))
colors_map = {
    "alpha_struct": "#2196F3",
    "alpha_spec": "#FF9800",
    "alpha_hyper": "#FF5722",
    "beta_struct": "#4CAF50",
    "beta_spec": "#9C27B0",
    "beta_hyper": "#E91E63",
}
bar_colors = [colors_map.get(summary_df.iloc[i]["category"].replace("Alpha~Struct", "alpha_struct")
              .replace("Alpha~Spec(RGB)", "alpha_spec").replace("Alpha~Spec(Hyper)", "alpha_hyper")
              .replace("Beta~Struct", "beta_struct").replace("Beta~Spec(RGB)", "beta_spec")
              .replace("Beta~Spec(Hyper)", "beta_hyper"), "#999")
              for i in range(len(summary_df))]
# Simpler approach
cat_to_color = {
    "Alpha~Struct": "#2196F3",
    "Alpha~Spec(RGB)": "#FF9800",
    "Alpha~Spec(Hyper)": "#FF5722",
    "Beta~Struct": "#4CAF50",
    "Beta~Spec(RGB)": "#9C27B0",
    "Beta~Spec(Hyper)": "#E91E63",
}
bar_colors = [cat_to_color[row["category"]] for _, row in summary_df.iterrows()]

bars = ax.bar(x, summary_df["sig_rate"] * 100, color=bar_colors, edgecolor="black", linewidth=0.5)

# Add mean r on top
for i, row in summary_df.iterrows():
    ax.text(i, row["sig_rate"] * 100 + 1,
            f"r={row['mean_r']:.2f}\nn={row['n_site_years']}",
            ha="center", va="bottom", fontsize=6.5)

ax.set_xticks(x)
ax.set_xticklabels([l.split("→")[0].strip() for l in x_labels], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Significant site-years (%)", fontsize=11)
ax.set_title("Summary: Within-Site Regression Significance Rate\n"
             "(% of site-years with p < 0.05)", fontsize=13, fontweight="bold")
ax.axhline(y=5, color="red", linestyle="--", alpha=0.5, label="5% chance level")

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, edgecolor="black", label=l)
                   for l, c in cat_to_color.items()]
ax.legend(handles=legend_elements, loc="upper right", fontsize=9)
ax.set_ylim(0, max(summary_df["sig_rate"] * 100) + 15)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
fig.savefig(OUT_DIR / "comprehensive_summary.png", dpi=150, bbox_inches="tight")
print(f"Saved: comprehensive_summary.png")

# ── Print final summary table ──────────────────────────────────────────────

print("\n" + "="*80)
print("COMPREHENSIVE SUMMARY TABLE")
print("="*80)
print(f"{'Label':<35s} {'N':>4s} {'mean r':>7s} {'mean R²':>7s} {'sig%':>5s} {'+sig':>4s} {'-sig':>4s}")
print("-"*80)
for cat_label in categories:
    cat_key = categories[cat_label]
    cat_rows = summary_df[summary_df["category"] == cat_label]
    if len(cat_rows) == 0:
        continue
    print(f"\n  [{cat_label}]")
    for _, row in cat_rows.iterrows():
        print(f"  {row['label']:<33s} {row['n_site_years']:4d} {row['mean_r']:7.3f} "
              f"{row['mean_r2']:7.3f} {row['sig_rate']*100:5.1f} {row['pos_sig']:4d} {row['neg_sig']:4d}")

# Save results
all_res_list = []
for label, (res, cat) in all_results.items():
    r = res.copy()
    r["label"] = label
    r["category"] = cat
    all_res_list.append(r)
combined = pd.concat(all_res_list, ignore_index=True)
combined.to_csv("E:/neon_lidar/model_results/yearly_regression_comprehensive.csv", index=False)
print(f"\nResults CSV: {len(combined)} rows saved")
print("\nDone.")
