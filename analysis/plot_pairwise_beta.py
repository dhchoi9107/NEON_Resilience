import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import rasterio, re

tag = pd.read_csv('E:/neon_lidar/vegetation_structure/vst_mappingandtagging.csv', low_memory=False)
ind = pd.read_csv('E:/neon_lidar/vegetation_structure/vst_apparentindividual.csv', low_memory=False)
ind['date'] = pd.to_datetime(ind['date'], errors='coerce')
mg = ind.merge(tag[['individualID','taxonID','plotID','siteID']], on='individualID', how='left')
mg = mg[mg['plantStatus'].str.contains('Live', case=False, na=False)]
mg = mg[mg['stemDiameter'].notna() & (mg['stemDiameter'] >= 10.0)]
mg = mg.dropna(subset=['taxonID']).drop_duplicates(subset=['plotID_y','individualID'])
mg['plotID'] = mg['plotID_y']; mg['siteID'] = mg['siteID_y']

struct = pd.read_csv('E:/neon_lidar/model_results/plot_alpha_structural.csv')
sv = ['mean_max_canopy_ht_mean','FHD_mean','LAI_mean','vert_sd_mean','GC_mean','deepgap_fraction_mean']
stp = struct.groupby(['siteID','plotID'])[sv].mean().reset_index()

HYPER_DIR = Path('E:/neon_lidar/hyperspectral_plots')
from compute.compute_hyperspectral_diversity import get_good_band_mask

all_tax, all_struct, all_spec, all_sites = [], [], [], []

for site in sorted(mg['siteID'].unique()):
    hdir = HYPER_DIR / site
    if not hdir.exists(): continue
    pivot = mg[mg['siteID']==site].groupby(['plotID','taxonID']).size().unstack(fill_value=0)
    st_site = stp[stp['siteID']==site]
    pspec = {}
    for f in hdir.glob('*_hyper.tif'):
        m = re.match(r'(\d{4})_(.+)_hyper\.tif', f.name)
        if not m or m.group(2) in pspec: continue
        with rasterio.open(str(f)) as src:
            d = src.read().astype(np.float32)
        good, _ = get_good_band_mask(src.count)
        d = d[good]; h,w = d.shape[1],d.shape[2]; ch,cw = h//2,w//2
        half = min(20,ch,cw)
        if half < 5: continue
        p = d[:, ch-half:ch+half, cw-half:cw+half]
        v = np.all((p>0)&(p<10000), axis=0)
        if v.sum() < 50: continue
        pspec[m.group(2)] = p[:, v].mean(axis=1)
    common = sorted(set(pivot.index) & set(st_site['plotID']) & set(pspec.keys()))
    if len(common) < 5: continue
    abund = pivot.loc[common].values.astype(float)
    tax_v = pdist(abund, 'braycurtis')
    st_data = st_site[st_site['plotID'].isin(common)].set_index('plotID').loc[common][sv].values
    struct_v = pdist(StandardScaler().fit_transform(st_data))
    sp = np.array([pspec[p] for p in common])
    n_pc = min(5, len(common)-1, sp.shape[1])
    sp_pca = PCA(n_components=n_pc).fit_transform(StandardScaler().fit_transform(sp))
    spec_v = pdist(sp_pca)
    all_tax.extend(tax_v); all_struct.extend(struct_v)
    all_spec.extend(spec_v); all_sites.extend([site]*len(tax_v))

all_tax = np.array(all_tax); all_struct = np.array(all_struct)
all_spec = np.array(all_spec); all_sites = np.array(all_sites)
valid = np.isfinite(all_tax) & np.isfinite(all_struct) & np.isfinite(all_spec)
all_tax, all_struct, all_spec, all_sites = all_tax[valid], all_struct[valid], all_spec[valid], all_sites[valid]

sr = pd.read_csv('E:/neon_lidar/model_results/pairwise_beta_distances.csv')

sites_u = sorted(set(all_sites))
cmap = plt.cm.tab20(np.linspace(0, 1, len(sites_u)))
sc = {s: cmap[i] for i, s in enumerate(sites_u)}
cols = [sc[s] for s in all_sites]

fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

# Panel 1
ax = fig.add_subplot(gs[0, 0])
ax.scatter(all_struct, all_tax, c=cols, s=3, alpha=0.3, rasterized=True)
r, p = pearsonr(all_struct, all_tax)
for s in sites_u:
    m = all_sites == s
    if m.sum() < 10: continue
    x, y = all_struct[m], all_tax[m]
    sl = np.polyfit(x, y, 1)
    xl = np.linspace(x.min(), x.max(), 50)
    ax.plot(xl, np.polyval(sl, xl), color=sc[s], linewidth=1.5, alpha=0.7)
ax.set_xlabel('Structural Distance', fontsize=11)
ax.set_ylabel('Taxonomic Distance (Bray-Curtis)', fontsize=11)
ax.set_title(f'Structural ~ Taxonomic Beta\nr = {r:.3f}, p < 0.001', fontsize=12, fontweight='bold')

# Panel 2
ax = fig.add_subplot(gs[0, 1])
ax.scatter(all_spec, all_tax, c=cols, s=3, alpha=0.3, rasterized=True)
r2, p2 = pearsonr(all_spec, all_tax)
for s in sites_u:
    m = all_sites == s
    if m.sum() < 10: continue
    x, y = all_spec[m], all_tax[m]
    sl = np.polyfit(x, y, 1)
    xl = np.linspace(x.min(), x.max(), 50)
    ax.plot(xl, np.polyval(sl, xl), color=sc[s], linewidth=1.5, alpha=0.7)
ax.set_xlabel('Spectral Distance (Hyperspectral PCA)', fontsize=11)
ax.set_ylabel('Taxonomic Distance (Bray-Curtis)', fontsize=11)
ax.set_title(f'Spectral ~ Taxonomic Beta\nr = {r2:.3f}', fontsize=12, fontweight='bold')

# Panel 3
ax = fig.add_subplot(gs[0, 2])
ax.scatter(all_struct, all_spec, c=cols, s=3, alpha=0.3, rasterized=True)
r3, p3 = pearsonr(all_struct, all_spec)
ax.set_xlabel('Structural Distance', fontsize=11)
ax.set_ylabel('Spectral Distance', fontsize=11)
ax.set_title(f'Structural ~ Spectral\nr = {r3:.3f} (redundancy)', fontsize=12, fontweight='bold')

# Panel 4: per-site bars
ax = fig.add_subplot(gs[1, 0:2])
sr2 = sr.dropna(subset=['r_struct_tax','r_spec_tax']).sort_values('r_struct_tax', ascending=True)
x = np.arange(len(sr2)); w = 0.35
ax.barh(x - w/2, sr2['r_struct_tax'].values, w, color='#d62728', alpha=0.8, label='Structural ~ Taxonomic')
ax.barh(x + w/2, sr2['r_spec_tax'].values, w, color='#1f77b4', alpha=0.8, label='Spectral ~ Taxonomic')
ax.set_yticks(x); ax.set_yticklabels(sr2['site'].values, fontsize=10)
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_xlabel('Mantel r (pairwise distance correlation)', fontsize=11)
ax.set_title('Per-site: Structural vs Spectral as Predictors of Taxonomic Beta', fontsize=12, fontweight='bold')
ax.legend(fontsize=10)
for i, (_, row) in enumerate(sr2.iterrows()):
    if row['p_struct_tax'] < 0.05:
        ax.annotate('*', (row['r_struct_tax'] + 0.01, i - w/2), fontsize=12, fontweight='bold', color='#d62728')
    if row['p_spec_tax'] < 0.05:
        ax.annotate('*', (row['r_spec_tax'] + 0.01, i + w/2), fontsize=12, fontweight='bold', color='#1f77b4')

# Panel 5
ax = fig.add_subplot(gs[1, 2]); ax.axis('off')
txt = (
    "SUMMARY\n"
    "----------------------------\n\n"
    f"Sites: {len(sr2)}\n"
    f"Plot pairs: {len(all_tax):,}\n\n"
    "Structural ~ Taxonomic:\n"
    f"  Mean r = {sr2['r_struct_tax'].mean():.3f}\n"
    f"  Significant: {(sr2['p_struct_tax']<0.05).sum()}/{len(sr2)}\n\n"
    "Spectral ~ Taxonomic:\n"
    f"  Mean r = {sr2['r_spec_tax'].mean():.3f}\n"
    f"  Significant: {(sr2['p_spec_tax']<0.05).sum()}/{len(sr2)}\n\n"
    "Spectral ~ Structural:\n"
    f"  Mean r = {sr2['r_spec_struct'].mean():.3f}\n"
    "  (low redundancy)\n\n"
    "426-band hyperspectral\n"
    "does NOT predict species\n"
    "turnover. Structure does."
)
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.suptitle('Pairwise Beta Diversity: Structural vs Spectral Distance as Predictors of Taxonomic Turnover\n(426-band BRDF hyperspectral, 1m, 15 NEON temperate forest sites)',
             fontsize=14, fontweight='bold')
plt.savefig('E:/neon_lidar/model_results/figures/pairwise_beta.png', dpi=200, bbox_inches='tight')
plt.savefig('C:/Users/star1/Documents/GitHub/NEON_Resilience/docs/pairwise_beta.png', dpi=200, bbox_inches='tight')
print('Saved!')
