"""
Framework figure (multivariate frame): (A) incremental R2 beyond domain across the nested
sequence State -> +Dynamics -> +Function, with LRT significance; (B) variance partition into
unique block fractions + shared. Reads O0_framework CSVs.
Out: papers/UNIFIED/figures/O0_framework/O0_framework.png
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
RES = os.path.join(BASE, "papers", "UNIFIED", "results", "O0_framework")
OUTF = os.path.join(BASE, "papers", "UNIFIED", "figures", "O0_framework"); os.makedirs(OUTF, exist_ok=True)

nm = pd.read_csv(os.path.join(RES, "nested_models.csv"))
vp = pd.read_csv(os.path.join(RES, "variance_partition.csv"))
cw = pd.read_csv(os.path.join(RES, "wildboot_blocks.csv")).set_index("response")  # wild cluster bootstrap p (few-cluster valid)
ORDER = ["Hill q1", "Hill q2", "LCBD turnover", "LCBD nestedness"]

def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "n.s."

# ---- assemble panel A (increments beyond domain) ----
rowsA = []
for r in ORDER:
    g = nm[nm.response == r]
    b1 = g[g.model == "M1_state"].R2_beyond_domain.values[0]
    b2 = g[g.model == "M2_+dynamics"].R2_beyond_domain.values[0]
    b3 = g[g.model == "M3_+function"].R2_beyond_domain.values[0]
    p_dyn = cw.loc[r, "p_dyn_wildboot"]      # wild cluster bootstrap (primary, few-cluster valid)
    p_fun = cw.loc[r, "p_fun_wildboot"]
    rowsA.append(dict(response=r, state=b1, dyn=b2-b1, fun=b3-b2, p_dyn=p_dyn, p_fun=p_fun))
A = pd.DataFrame(rowsA)

fig, ax = plt.subplots(1, 2, figsize=(15, 5.6), constrained_layout=True)
x = np.arange(len(ORDER)); C = {"state":"#1f78b4", "dyn":"#33a02c", "fun":"#ff7f00"}
# Panel A: stacked incremental R2 beyond domain
ax[0].bar(x, A.state, color=C["state"], label="State (structure + spectral)")
ax[0].bar(x, A.dyn, bottom=A.state, color=C["dyn"], label="+ Dynamics (LiDAR trends)")
ax[0].bar(x, A.fun, bottom=A.state + A.dyn, color=C["fun"], label="+ Function (GPP)")
for i, rr in A.iterrows():
    ax[0].text(i, rr.state/2, f"{rr.state:.3f}", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax[0].text(i, rr.state+rr.dyn/2, f"+{rr.dyn:.3f}\n{star(rr.p_dyn)}", ha="center", va="center", color="white", fontsize=7.5)
    ytop = rr.state+rr.dyn+rr.fun
    ax[0].text(i, ytop+0.004, f"+{rr.fun:.3f} {star(rr.p_fun)}", ha="center", va="bottom", color=C["fun"], fontsize=7.5)
ax[0].set_xticks(x); ax[0].set_xticklabels(ORDER, fontsize=10)
ax[0].set_ylabel("Incremental R² beyond domain", fontsize=11)
ax[0].set_title("(A) Each dimension adds information sequentially", fontsize=12, fontweight="bold")
ax[0].legend(fontsize=9, loc="upper right"); ax[0].set_ylim(0, max(A.state+A.dyn+A.fun)*1.25)

# Panel B: variance partition (unique + shared)
blocks = ["unique_structure", "unique_spectral", "unique_dynamics", "unique_productivity", "shared_RS"]
labs = ["structure (unique)", "spectral (unique)", "dynamics (unique)", "productivity (unique)", "shared"]
cols = ["#1f78b4", "#a6cee3", "#33a02c", "#ff7f00", "#bdbdbd"]
vpO = vp.set_index("response").loc[ORDER]
bottom = np.zeros(len(ORDER))
for b, lab, c in zip(blocks, labs, cols):
    vals = vpO[b].values
    ax[1].bar(x, vals, bottom=bottom, color=c, label=lab); bottom += vals
ax[1].set_xticks(x); ax[1].set_xticklabels(ORDER, fontsize=10)
ax[1].set_ylabel("R² beyond domain (partitioned)", fontsize=11)
ax[1].set_title("(B) Structure carries the largest unique share;\nspectral ≈ 0 for beta", fontsize=12, fontweight="bold")
ax[1].legend(fontsize=8.5, loc="upper right")
_na = int(nm[nm.response == "Hill q1"]["n"].iloc[0]); _nb = int(nm[nm.response == "LCBD turnover"]["n"].iloc[0])
_ns = pd.read_csv(os.path.join(BASE, "data", "FINAL_v2_pooled_26.csv")).siteID.nunique()
fig.suptitle(f"Multivariate framework: complementary remote-sensing dimensions of tree diversity (n={_na} alpha, {_nb} beta; OLS + domain; stars = wild cluster bootstrap, {_ns} sites)", fontsize=11.5, fontweight="bold")
out = os.path.join(OUTF, "O0_framework.png")
fig.savefig(out, dpi=200, bbox_inches="tight"); print("saved", out)
print(A.round(4).to_string(index=False))
