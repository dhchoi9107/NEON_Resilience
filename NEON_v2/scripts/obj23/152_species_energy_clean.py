"""
O3 clean figure: species-energy relationship. Tree diversity (Hill q1) vs independent GPP
(MODIS MOD17, PML-V2, NEON eddy-covariance towers). Monotonic positive; NO DHI / no hump framing.
Reads pre-extracted GPP (no GEE needed).
Out: papers/UNIFIED/figures/O3_productivity/O3_species_energy.png
     papers/UNIFIED/results/O3_productivity/species_energy_stats.csv
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D = os.path.join(BASE, "data")
OUTF = os.path.join(BASE, "papers", "UNIFIED", "figures", "O3_productivity")
OUTR = os.path.join(BASE, "papers", "UNIFIED", "results", "O3_productivity")
os.makedirs(OUTF, exist_ok=True); os.makedirs(OUTR, exist_ok=True)

# --- load diversity + GPP proxies ---
df = pd.read_csv(os.path.join(D, "FINAL_v2_pooled.csv"))
df = df.merge(pd.read_csv(os.path.join(D, "plot_modis_gpp.csv")), on="plotID", how="left")
df = df.merge(pd.read_csv(os.path.join(D, "plot_pml_gpp.csv")), on="plotID", how="left")
d = df[df.sample_coverage >= 0.9].copy()

tower = pd.read_csv(os.path.join(D, "site_tower_gpp.csv"))  # siteID, tower_gpp (annual gC/m2/yr)

def zfit(x, y):
    """standardized linear + quadratic test; return slope beta, r, p_lin, quad beta/p."""
    m = pd.DataFrame({"x": x, "y": y}).dropna()
    zx = (m.x - m.x.mean()) / m.x.std()
    lin = smf.ols("y ~ zx", pd.DataFrame({"y": m.y, "zx": zx})).fit()
    quad = smf.ols("y ~ zx + I(zx**2)", pd.DataFrame({"y": m.y, "zx": zx})).fit()
    qk = [c for c in quad.params.index if "** 2" in c][0]
    r, p = st.pearsonr(m.x, m.y)
    return dict(n=len(m), beta=lin.params["zx"], p_lin=lin.pvalues["zx"],
                r=r, p_r=p, quad_beta=quad.params[qk], quad_p=quad.pvalues[qk],
                xmin=m.x.min(), xmax=m.x.max())

panels = [
    ("modis_gpp", "MODIS MOD17 GPP", "plot", d, "Hill_q1"),
    ("pml_gpp",   "PML-V2 GPP",       "plot", d, "Hill_q1"),
    ("tower_gpp", "NEON tower GPP",   "site", None, "Hill_q1"),
]
# site-mean diversity for tower panel
site_div = d.groupby("siteID").agg(Hill_q1=("Hill_q1", "mean")).reset_index()
tw = tower.merge(site_div, on="siteID", how="inner")

GREEN = "#00695c"; RED = "#c62828"; UNITS = {"modis_gpp": "gC m$^{-2}$ d$^{-1}$",
         "pml_gpp": "gC m$^{-2}$ d$^{-1}$", "tower_gpp": "gC m$^{-2}$ yr$^{-1}$"}
rows = []
fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
for ax, (col, lab, scale, data, yv) in zip(axes, panels):
    if scale == "plot":
        dd = data[[col, yv, "siteID"]].dropna()
        ax.scatter(dd[col], dd[yv], s=12, alpha=0.25, color=GREEN, edgecolor="none",
                   label="plot (n=%d)" % len(dd))
        sm = dd.groupby("siteID").agg(x=(col, "median"), y=(yv, "mean"))
        ax.scatter(sm.x, sm.y, s=55, color=RED, zorder=3, edgecolor="white",
                   linewidth=0.6, label="site mean")
        s = zfit(dd[col], dd[yv]); s_site = zfit(sm.x, sm.y)
        rlabel = f"plot r={s['r']:+.2f}***\nsite r={s_site['r']:+.2f}" if s['p_r'] < .001 \
                 else f"plot r={s['r']:+.2f} (p={s['p_r']:.2g})\nsite r={s_site['r']:+.2f}"
        xr = np.linspace(dd[col].min(), dd[col].max(), 50)
    else:
        ax.scatter(tw[col], tw[yv], s=90, color=RED, zorder=3, edgecolor="white",
                   linewidth=0.8, label="site (n=%d)" % len(tw))
        for _, r_ in tw.iterrows():
            ax.annotate(r_.siteID, (r_[col], r_[yv]), fontsize=6, alpha=0.6,
                        xytext=(3, 3), textcoords="offset points")
        s = zfit(tw[col], tw[yv])
        rlabel = f"r={s['r']:+.2f}" + ("**" if s['p_r'] < .01 else "*" if s['p_r'] < .05
                 else f" (p={s['p_r']:.2g})")
        xr = np.linspace(tw[col].min(), tw[col].max(), 50)
        dd = tw.rename(columns={col: col}); sm = None
    # linear fit line (raw units)
    src = dd if scale == "plot" else tw
    b1, b0 = np.polyfit(src[col].dropna(), src.loc[src[col].notna(), yv], 1)
    ax.plot(xr, b0 + b1 * xr, color="black", lw=2, zorder=4)
    ax.set_xlabel(f"GPP ({UNITS[col]})", fontsize=11)
    ax.set_title(lab, fontsize=12, fontweight="bold")
    is_hump = (s["quad_beta"] < 0) and (s["quad_p"] < .05)   # inverted-U only
    quadtxt = "hump (inverted-U)" if is_hump else "monotonic, no hump"
    ax.text(0.04, 0.96, rlabel + f"\nslope>0; {quadtxt}", transform=ax.transAxes,
            va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.legend(loc="lower right", fontsize=8, framealpha=0.8)
    rows.append(dict(proxy=lab, scale=scale, **s))

axes[0].set_ylabel("Tree diversity (Hill q1)", fontsize=11)
fig.suptitle("Species–energy: tree diversity increases monotonically with independent GPP",
             fontsize=13, fontweight="bold")
out = os.path.join(OUTF, "O3_species_energy.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); print("saved", out)
pd.DataFrame(rows).to_csv(os.path.join(OUTR, "species_energy_stats.csv"), index=False)
print(pd.DataFrame(rows)[["proxy", "scale", "n", "beta", "p_lin", "r", "p_r", "quad_beta", "quad_p"]].to_string(index=False))
