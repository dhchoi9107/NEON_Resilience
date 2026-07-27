"""
Build pooled predictor rows for the 7 NEW sites, matching lidar_pooled_predictors.csv
schema EXACTLY, then concatenate with the existing 19-site file -> 26-site version.

Replicates the ORIGINAL methods verbatim:
  - alpha: sampling_completeness.hill_coverage (Chao-Jost coverage, Hill q1/q2, rare10)
  - beta (non-rare): beta_coverage_filtered (within-site Baselga, cov>=0.9)
  - beta (rare): rarefied_lcbd (within-site Baselga, cov>=0.9 & N>=10, RARE_N=10, 30 draws, seed 11)
  - structure: per-plot mean/sd/trend/nyears over per_year_newsites_10m years
LCBD is WITHIN-SITE => existing 19 sites are unchanged; we only add the 7 new sites.

Output (NEW file, never overwrites): NEON_v2/data/lidar_pooled_predictors_26.csv
"""
import os, numpy as np, pandas as pd
from math import comb
from scipy.spatial.distance import pdist, squareform

ROOT = r"C:\Users\star1\Documents\GitHub\NEON_Resilience"
D = os.path.join(ROOT, "NEON_v2", "data")
VST = "E:/neon_lidar/vegetation_structure"
NEW = ["DELA", "LENO", "UKFS", "YELL", "BONA", "DEJU", "HEAL"]
DOM = {"DELA": "D08", "LENO": "D08", "UKFS": "D06", "YELL": "D12",
       "BONA": "D19", "DEJU": "D19", "HEAL": "D19"}
STRUCT = ['Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV',
          'Gini','VCI','FHD','LAI','Q95','Ht_Ratio']
RARE_N = 10; N_DRAWS = 30; RNG = np.random.RandomState(11)

# ── 1. pooled SAD (DBH>=10, live, valid taxon), identical filters ──
mt = pd.read_csv(f"{VST}/vst_mappingandtagging.csv", usecols=['individualID','plotID','taxonID','siteID'], low_memory=False)
ai = pd.read_csv(f"{VST}/vst_apparentindividual.csv", usecols=['individualID','plantStatus','stemDiameter','siteID'], low_memory=False)
df = ai.merge(mt[['individualID','plotID','taxonID']], on='individualID', how='left')
df = df[df['plantStatus'].astype(str).str.contains('Live', na=False)]
df = df[df['stemDiameter'] >= 10]
df = df[df['taxonID'].notna() & ~df['taxonID'].astype(str).str.contains('2PLANT|UNK', na=False)]
df = df.drop_duplicates('individualID')
df = df[df['siteID'].isin(NEW)]
print(f"new-site individuals: {len(df)}, plots: {df['plotID'].nunique()}", flush=True)

# ── 2. alpha: Chao-Jost coverage + Hill q1/q2 (sampling_completeness.hill_coverage) ──
def hill_coverage(counts):
    counts = np.asarray([c for c in counts if c > 0]); N = counts.sum(); S = len(counts)
    if N == 0: return dict(S=0, N=0, q1=0, q2=0, cov=np.nan)
    p = counts / N
    q1 = float(np.exp(-(p*np.log(p)).sum()))
    q2 = float(1.0/(p**2).sum())
    f1 = int((counts == 1).sum()); f2 = int((counts == 2).sum())
    cov = 1 - (f1/N)*(((N-1)*f1)/((N-1)*f1 + 2*max(f2,1))) if N > 1 else np.nan
    return dict(S=S, N=int(N), q1=q1, q2=q2, cov=cov)

arows = []
for plot, g in df.groupby('plotID'):
    counts = g.groupby('taxonID').size().values
    arows.append(dict(plotID=plot, siteID=g['siteID'].iloc[0], **hill_coverage(counts)))
alpha = pd.DataFrame(arows)
print(f"alpha plots: {len(alpha)} | cov>=0.9: {(alpha['cov']>=0.9).sum()}", flush=True)

# ── 3. beta non-rare (within-site Baselga, cov>=0.9) ──
def baselga(pa):
    n = pa.shape[0]; sim = np.zeros((n,n)); sne = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j: continue
            a = np.sum((pa[i]>0)&(pa[j]>0)); b = np.sum((pa[i]>0)&(pa[j]==0)); c = np.sum((pa[i]==0)&(pa[j]>0))
            sor = (b+c)/(2*a+b+c) if (2*a+b+c) > 0 else 0
            s = min(b,c)/(a+min(b,c)) if (a+min(b,c)) > 0 else 0
            sim[i,j] = s; sne[i,j] = sor - s
    return sim, sne

good = set(alpha[alpha['cov'] >= 0.9]['plotID'])
brows = []
for site, g in df.groupby('siteID'):
    g = g[g['plotID'].isin(good)]
    piv = g.groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    if piv.shape[0] < 2: continue
    ab = piv.values.astype(float); pa = (ab > 0).astype(float); plots = piv.index.tolist(); n = len(plots)
    BC = squareform(pdist(ab, metric='braycurtis')); sim, sne = baselga(pa)
    for i, p in enumerate(plots):
        oth = np.arange(n) != i
        brows.append(dict(plotID=p, LCBD_bray=BC[i][oth].mean(),
                          LCBD_turnover=sim[i][oth].mean(), LCBD_nestedness=sne[i][oth].mean()))
beta = pd.DataFrame(brows)

# ── 4. beta rarefied (within-site, cov>=0.9 & N>=10, RARE_N=10, 30 draws) ──
def rarefy_vec(counts, m):
    counts = np.asarray(counts).astype(int); N = counts.sum()
    if N <= m: return counts.astype(float)
    pool = np.repeat(np.arange(len(counts)), counts)
    pick = RNG.choice(pool, m, replace=False)
    return np.bincount(pick, minlength=len(counts)).astype(float)

goodr = set(alpha[(alpha['cov'] >= 0.9) & (alpha['N'] >= RARE_N)]['plotID'])
rrows = []
for site, g in df.groupby('siteID'):
    g = g[g['plotID'].isin(goodr)]
    piv = g.groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    if piv.shape[0] < 2: continue
    plots = piv.index.tolist(); base = piv.values.astype(float); n = len(plots)
    ab_ = np.zeros(n); tu = np.zeros(n); ne = np.zeros(n)
    for _ in range(N_DRAWS):
        rar = np.vstack([rarefy_vec(base[i], RARE_N) for i in range(n)])
        BC = squareform(pdist(rar, metric='braycurtis')); pa = (rar > 0).astype(float); sim, sne = baselga(pa)
        for i in range(n):
            oth = np.arange(n) != i
            ab_[i] += BC[i][oth].mean(); tu[i] += sim[i][oth].mean(); ne[i] += sne[i][oth].mean()
    for i, p in enumerate(plots):
        rrows.append(dict(plotID=p, LCBD_bray_rare=ab_[i]/N_DRAWS,
                          LCBD_turnover_rare=tu[i]/N_DRAWS, LCBD_nestedness_rare=ne[i]/N_DRAWS))
betar = pd.DataFrame(rrows)
print(f"beta non-rare plots: {len(beta)} | rarefied plots: {len(betar)}", flush=True)

# ── 5. structure pooled (mean/sd/trend/nyears) from per_year_newsites_10m ──
py = pd.read_csv("scripts_pipeline/_pipeline_state/per_year_newsites_10m.csv")
def slope(y, x):
    ok = np.isfinite(y)
    return float(np.polyfit(x[ok], y[ok], 1)[0]) if ok.sum() >= 3 else np.nan
srows = []
for pid, g in py.groupby('plotID'):
    g = g.sort_values('year'); yrs = g['year'].values.astype(float)
    rec = {'plotID': pid, 'siteID': g['siteID'].iloc[0]}
    for m in STRUCT:
        v = g[m].values.astype(float); ok = np.isfinite(v)
        rec[f"{m}_mean"] = float(np.nanmean(v)) if ok.sum() else np.nan
        rec[f"{m}_sd"] = float(np.nanstd(v)) if ok.sum() >= 2 else np.nan
        rec[f"{m}_trend"] = slope(v, yrs) if ok.sum() >= 3 else np.nan
        rec[f"{m}_nyears"] = int(ok.sum())
    srows.append(rec)
struct = pd.DataFrame(srows)
print(f"structure plots: {len(struct)}", flush=True)

# ── 6. assemble new-site rows to match lidar_pooled_predictors schema ──
new = (struct
       .merge(alpha.rename(columns={'q1':'Hill_q1','q2':'Hill_q2','cov':'sample_coverage',
                                    'S':'richness_pooled','N':'abundance_pooled'})[
           ['plotID','Hill_q1','Hill_q2','sample_coverage','richness_pooled','abundance_pooled']], on='plotID', how='left')
       .merge(beta, on='plotID', how='left')
       .merge(betar, on='plotID', how='left'))
new['domain'] = new['siteID'].map(DOM)

# align columns to existing file, concat
existing = pd.read_csv(os.path.join(D, "lidar_pooled_predictors.csv"))
new = new.reindex(columns=existing.columns)
comb = pd.concat([existing, new], ignore_index=True)
OUT = os.path.join(D, "lidar_pooled_predictors_26.csv")
comb.to_csv(OUT, index=False)
print(f"\n-> {OUT}", flush=True)
print(f"total: {len(comb)} plots, {comb.siteID.nunique()} sites", flush=True)
print(f"new-site rows: {len(new)} | Hill_q1 nn: {new.Hill_q1.notna().sum()} | "
      f"cov>=0.9: {(new.sample_coverage>=0.9).sum()} | LCBD_turn_rare nn: {new.LCBD_turnover_rare.notna().sum()}", flush=True)
print("new-site coverage>=0.9 plots by site:", flush=True)
print(new[new.sample_coverage>=0.9].groupby('siteID').size().to_string(), flush=True)
print("\nsanity (new-site medians):", flush=True)
print(new.groupby('siteID').agg(Ht=('Canopy_Ht_mean','median'), Hill1=('Hill_q1','median'),
      turn=('LCBD_turnover_rare','median')).round(2).to_string(), flush=True)
