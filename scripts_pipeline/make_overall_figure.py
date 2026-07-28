"""
Single overall-results figure (26 sites + TREE recovered) summarizing the whole story:
 (A) sequential incremental R2 beyond domain (State -> +Dynamics -> +Function), wildboot stars
 (B) variance partition into unique block contributions + shared
 (C) fully-standardized coefficient heatmap (predictors x 4 diversity responses)
 (D) context modulation: GPP-alpha slope steepens with stand age; structure-alpha slope
     weakens with disturbance severity (simple slopes at -1/0/+1 SD)
Reads O0_framework CSVs (current). Out: papers/UNIFIED/figures/OVERALL_results.png
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
RES = os.path.join(BASE, "papers", "UNIFIED", "results", "O0_framework")
OUT = os.path.join(BASE, "papers", "UNIFIED", "figures", "OVERALL_results.png")

nm = pd.read_csv(os.path.join(RES, "nested_models.csv"))
vp = pd.read_csv(os.path.join(RES, "variance_partition.csv")).set_index("response")
wb = pd.read_csv(os.path.join(RES, "wildboot_blocks.csv")).set_index("response")
cf = pd.read_csv(os.path.join(RES, "coeffs_fullstd.csv"))
ss = pd.read_csv(os.path.join(RES, "simple_slopes.csv"))
ns = pd.read_csv(os.path.join(BASE, "data", "FINAL_v2_pooled_26.csv")).siteID.nunique()
na = int(nm[nm.response == "Hill q1"]["n"].iloc[0]); nb = int(nm[nm.response == "LCBD turnover"]["n"].iloc[0])

RESP = ["Hill q1", "Hill q2", "LCBD turnover", "LCBD nestedness"]
LAB = ["Hill q1\n(alpha)", "Hill q2\n(alpha)", "LCBD\nturnover", "LCBD\nnestedness"]
def star(p): return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "n.s."

def r2beyond(resp, model):
    return float(nm[(nm.response == resp) & (nm.model == model)]["R2_beyond_domain"].iloc[0])

fig, axes = plt.subplots(2, 2, figsize=(15, 11))
fig.suptitle(f"Remote-sensing dimensions of tree diversity across NEON — overall results "
             f"(n={na} alpha / {nb} beta plots, {ns} sites; stars = wild cluster bootstrap)",
             fontsize=13.5, fontweight="bold")

# ---------- (A) sequential incremental R2 ----------
ax = axes[0, 0]; x = np.arange(len(RESP)); w = 0.6
state = [r2beyond(r, "M1_state") for r in RESP]
dyn = [r2beyond(r, "M2_+dynamics") - r2beyond(r, "M1_state") for r in RESP]
fun = [r2beyond(r, "M3_+function") - r2beyond(r, "M2_+dynamics") for r in RESP]
ax.bar(x, state, w, label="State (structure + spectral)", color="#2c7fb8")
ax.bar(x, dyn, w, bottom=state, label="+ Dynamics (LiDAR trends)", color="#41ab5d")
ax.bar(x, fun, w, bottom=np.array(state) + np.array(dyn), label="+ Function (GPP)", color="#f16913")
for i, r in enumerate(RESP):
    ax.text(i, state[i] / 2, f"{state[i]:.3f}", ha="center", va="center", color="white", fontweight="bold", fontsize=9)
    ax.text(i, state[i] + dyn[i] / 2, f"+{dyn[i]:.3f}\n{star(wb.loc[r,'p_dyn_wildboot'])}", ha="center", va="center", color="white", fontsize=8)
    top = state[i] + dyn[i] + fun[i]
    ax.text(i, top + 0.004, f"+{fun[i]:.3f} {star(wb.loc[r,'p_fun_wildboot'])}", ha="center", va="bottom", color="#f16913", fontsize=8.5)
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylabel("Incremental R² beyond domain")
ax.set_title("(A) Each dimension adds information sequentially", fontweight="bold", fontsize=11)
ax.legend(fontsize=8.5, loc="upper right"); ax.set_ylim(0, 0.20)

# ---------- (B) variance partition ----------
ax = axes[0, 1]
blocks = [("unique_structure", "structure", "#08519c"), ("unique_spectral", "spectral", "#9ecae1"),
          ("unique_dynamics", "dynamics", "#41ab5d"), ("unique_productivity", "productivity", "#f16913"),
          ("shared_RS", "shared", "#bdbdbd")]
bottom = np.zeros(len(RESP))
for col, lab, c in blocks:
    vals = [vp.loc[r, col] for r in RESP]
    ax.bar(x, vals, w, bottom=bottom, label=lab, color=c); bottom += np.array(vals)
ax.set_xticks(x); ax.set_xticklabels(LAB); ax.set_ylabel("R² beyond domain (partitioned)")
ax.set_title("(B) Unique contribution of each dimension\n(structure→alpha; dynamics→nestedness; turnover carries little unique signal)", fontweight="bold", fontsize=11)
ax.legend(fontsize=8.5, loc="upper right")

# ---------- (C) coefficient heatmap ----------
ax = axes[1, 0]
PRED = ["VCI_mean", "LAI_mean", "Rugosity_mean", "Vert_CV_mean", "EVI_mean",
        "Rumple_trend", "FHD_trend", "LAI_trend", "Ht_Ratio_trend", "VCI_trend",
        "modis_gpp", "pml_gpp"]
PLAB = ["VCI (struct)", "LAI (struct)", "Rugosity (struct)", "Vert-CV (struct)", "EVI (spectral)",
        "Rumple trend (dyn)", "FHD trend (dyn)", "LAI trend (dyn)", "Ht-ratio trend (dyn)", "VCI trend (dyn)",
        "MODIS GPP (fun)", "PML GPP (fun)"]
M = np.full((len(PRED), len(RESP)), np.nan); P = np.full_like(M, np.nan)
for j, r in enumerate(RESP):
    sub = cf[cf.response == r].set_index("predictor")
    for i, p in enumerate(PRED):
        if p in sub.index:
            M[i, j] = sub.loc[p, "beta_fullstd"]; P[i, j] = sub.loc[p, "p_clustered"]
norm = TwoSlopeNorm(vmin=-0.4, vcenter=0, vmax=0.4)
im = ax.imshow(M, cmap="RdBu_r", norm=norm, aspect="auto")
ax.set_xticks(range(len(RESP))); ax.set_xticklabels(LAB, fontsize=9)
ax.set_yticks(range(len(PRED))); ax.set_yticklabels(PLAB, fontsize=9)
for i in range(len(PRED)):
    for j in range(len(RESP)):
        if not np.isnan(M[i, j]):
            s = "***" if P[i,j] < 1e-3 else "**" if P[i,j] < 1e-2 else "*" if P[i,j] < 5e-2 else ""
            ax.text(j, i, f"{M[i,j]:+.2f}{s}", ha="center", va="center", fontsize=7.5,
                    color="white" if abs(M[i,j]) > 0.22 else "black")
ax.set_title("(C) Standardized coefficients (full model)\n* p<.05  ** p<.01  *** p<.001 (site-clustered)", fontweight="bold", fontsize=11)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="β (standardized)")

# ---------- (D) context modulation ----------
ax = axes[1, 1]
z = [-1, 0, 1]
g = ss[(ss.rp == "Hill q1") & (ss.rs == "modis_gpp")].sort_values("ctx_z")
d = ss[(ss.rp == "Hill q1") & (ss.rs == "VCI_mean")].sort_values("ctx_z")
ax.errorbar(z, g.slope, yerr=g.se, marker="o", capsize=4, color="#f16913", lw=2, label="GPP→alpha  vs  stand age")
ax.errorbar(z, d.slope, yerr=d.se, marker="s", capsize=4, color="#08519c", lw=2, label="Structure(VCI)→alpha  vs  disturbance")
ax.axhline(0, color="grey", lw=0.8, ls="--")
for zz, sl, pv in zip(z, g.slope, g.p): ax.annotate(("*" if pv<0.05 else ""), (zz, sl), textcoords="offset points", xytext=(6,4), color="#f16913", fontsize=12)
for zz, sl, pv in zip(z, d.slope, d.p): ax.annotate(("*" if pv<0.05 else ""), (zz, sl), textcoords="offset points", xytext=(6,4), color="#08519c", fontsize=12)
ax.set_xticks(z); ax.set_xticklabels(["−1 SD\n(young / undisturbed)", "mean", "+1 SD\n(old / disturbed)"], fontsize=9)
ax.set_ylabel("simple slope on alpha diversity"); ax.set_xlabel("moderator (context)")
ax.set_title("(D) Context modulates the associations\nGPP–age trend suggestive (interaction n.s. at 26 sites); structure effect weakens with disturbance", fontweight="bold", fontsize=10)
ax.legend(fontsize=8.5, loc="upper left")

fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("saved", OUT)
