"""
BRDF-corrected Hyperspectral Analysis
========================================
1) BRDF corrected vs uncorrected spectral diversity comparison
2) Time-series pairwise beta (spectral vs taxonomic distance)
3) Disturbance effects on spectral change

Uses .001 BRDF-corrected (2013-2021) + .002 original BRDF (2022-2025)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
import rasterio
import re
import warnings
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from compute.compute_hyperspectral_diversity import get_good_band_mask

warnings.filterwarnings("ignore")

BRDF_DIR = Path("E:/neon_lidar/hyperspectral_brdf_corrected")
NOBRDF_DIR = Path("E:/neon_lidar/hyperspectral_plots")  # .002 (already BRDF by NEON)
NOBRDF_001_DIR = Path("E:/neon_lidar/hyperspectral_plots_001")  # .001 uncorrected
FIG_DIR = Path("E:/neon_lidar/model_results/figures")
OUT_DIR = Path("E:/neon_lidar/model_results")


def load_plot_mean_spectrum(tif_path, plot_half=20):
    """Load hyperspectral clip, return mean spectrum for center 40x40m."""
    with rasterio.open(str(tif_path)) as src:
        data = src.read().astype(np.float32)
        n_bands = src.count
    good, wl = get_good_band_mask(n_bands)
    data = data[good]
    h, w = data.shape[1], data.shape[2]
    ch, cw = h // 2, w // 2
    half = min(plot_half, ch, cw)
    if half < 5:
        return None
    patch = data[:, ch-half:ch+half, cw-half:cw+half]
    valid = np.all((patch > 0) & (patch < 10000), axis=0)
    if valid.sum() < 50:
        return None
    return patch[:, valid].mean(axis=1)


def discover_clips(base_dir, suffix="_hyper_brdf.tif"):
    """Discover all plot clips: {(site, plotID, year): path}"""
    clips = {}
    for site_dir in sorted(base_dir.iterdir()):
        if not site_dir.is_dir():
            continue
        site = site_dir.name
        for f in site_dir.glob(f"*{suffix}"):
            m = re.match(r"(\d{4})_(.+?)(?:_hyper(?:_brdf)?)?\.tif", f.name)
            if not m:
                continue
            year, plot_id = int(m.group(1)), m.group(2)
            clips[(site, plot_id, year)] = f
    return clips


# ═══════════════════════════════════════════════════════════════
# 1) BRDF corrected vs uncorrected comparison
# ═══════════════════════════════════════════════════════════════

def analysis_1_brdf_comparison():
    """Compare spectral diversity before/after BRDF correction."""
    print("=== Analysis 1: BRDF correction effect ===")

    brdf_clips = discover_clips(BRDF_DIR, "_hyper_brdf.tif")
    nobrdf_clips = discover_clips(NOBRDF_001_DIR, "_hyper.tif")

    # Find common (site, plot, year)
    common = sorted(set(brdf_clips.keys()) & set(nobrdf_clips.keys()))
    print(f"  Common plot-years: {len(common)}")

    if len(common) < 20:
        print("  Not enough common data for comparison")
        return

    # Per-site: compare pairwise spectral distances
    results = []
    for site in sorted(set(k[0] for k in common)):
        site_common = [(s, p, y) for s, p, y in common if s == site]
        if len(site_common) < 5:
            continue

        # Pick one year with most plots
        year_counts = {}
        for s, p, y in site_common:
            year_counts[y] = year_counts.get(y, 0) + 1
        best_year = max(year_counts, key=year_counts.get)
        plots = [(s, p, y) for s, p, y in site_common if y == best_year]

        if len(plots) < 5:
            continue

        # Load spectra
        brdf_spectra, nobrdf_spectra, plot_ids = [], [], []
        for s, p, y in plots:
            bs = load_plot_mean_spectrum(brdf_clips[(s, p, y)])
            ns = load_plot_mean_spectrum(nobrdf_clips[(s, p, y)])
            if bs is not None and ns is not None:
                brdf_spectra.append(bs)
                nobrdf_spectra.append(ns)
                plot_ids.append(p)

        if len(brdf_spectra) < 5:
            continue

        brdf_arr = np.array(brdf_spectra)
        nobrdf_arr = np.array(nobrdf_spectra)

        # Pairwise spectral distance
        brdf_pca = PCA(n_components=min(5, len(brdf_arr)-1)).fit_transform(
            StandardScaler().fit_transform(brdf_arr))
        nobrdf_pca = PCA(n_components=min(5, len(nobrdf_arr)-1)).fit_transform(
            StandardScaler().fit_transform(nobrdf_arr))

        brdf_dist = pdist(brdf_pca)
        nobrdf_dist = pdist(nobrdf_pca)

        r_dist, p_dist = pearsonr(brdf_dist, nobrdf_dist)

        # CV of spectra (within-site spectral variation)
        brdf_cv = np.mean(np.std(brdf_arr, axis=0) / (np.mean(brdf_arr, axis=0) + 1e-6))
        nobrdf_cv = np.mean(np.std(nobrdf_arr, axis=0) / (np.mean(nobrdf_arr, axis=0) + 1e-6))

        results.append({
            'site': site, 'year': best_year, 'n_plots': len(brdf_spectra),
            'r_distances': r_dist, 'brdf_cv': brdf_cv, 'nobrdf_cv': nobrdf_cv,
            'cv_ratio': brdf_cv / nobrdf_cv if nobrdf_cv > 0 else np.nan,
        })
        print(f"  {site} ({best_year}): {len(brdf_spectra)} plots, "
              f"dist r={r_dist:.3f}, CV ratio={brdf_cv/nobrdf_cv:.3f}")

    if results:
        rdf = pd.DataFrame(results)
        print(f"\n  Mean distance correlation (BRDF vs no-BRDF): {rdf['r_distances'].mean():.3f}")
        print(f"  Mean CV ratio (BRDF/no-BRDF): {rdf['cv_ratio'].mean():.3f}")
        rdf.to_csv(OUT_DIR / "brdf_comparison.csv", index=False)


# ═══════════════════════════════════════════════════════════════
# 2) Time-series pairwise beta
# ═══════════════════════════════════════════════════════════════

def analysis_2_timeseries_beta():
    """Pairwise spectral vs taxonomic beta across years."""
    print("\n=== Analysis 2: Time-series pairwise beta ===")

    brdf_clips = discover_clips(BRDF_DIR, "_hyper_brdf.tif")

    # Also add .002 clips (2022+)
    neon_brdf_clips = discover_clips(NOBRDF_DIR, "_hyper.tif")
    all_clips = {**brdf_clips, **neon_brdf_clips}  # .002 overwrites if overlap

    # Load taxonomy
    tag = pd.read_csv('E:/neon_lidar/vegetation_structure/vst_mappingandtagging.csv', low_memory=False)
    ind = pd.read_csv('E:/neon_lidar/vegetation_structure/vst_apparentindividual.csv', low_memory=False)
    ind['date'] = pd.to_datetime(ind['date'], errors='coerce')
    mg = ind.merge(tag[['individualID','taxonID','plotID','siteID']], on='individualID', how='left')
    mg = mg[mg['plantStatus'].str.contains('Live', case=False, na=False)]
    mg = mg[mg['stemDiameter'].notna() & (mg['stemDiameter'] >= 10.0)]
    mg = mg.dropna(subset=['taxonID']).drop_duplicates(subset=['plotID_y','individualID'])
    mg['plotID'] = mg['plotID_y']; mg['siteID'] = mg['siteID_y']

    # Structural data
    struct = pd.read_csv('E:/neon_lidar/model_results/plot_alpha_structural.csv')
    sv = ['mean_max_canopy_ht_mean','FHD_mean','LAI_mean','vert_sd_mean','GC_mean','deepgap_fraction_mean']
    stp = struct.groupby(['siteID','plotID'])[sv].mean().reset_index()

    results = []
    for site in sorted(set(k[0] for k in all_clips.keys())):
        site_keys = [(s, p, y) for s, p, y in all_clips if s == site]
        years = sorted(set(y for _, _, y in site_keys))

        # Taxonomy (pooled)
        site_sp = mg[mg['siteID'] == site]
        pivot = site_sp.groupby(['plotID', 'taxonID']).size().unstack(fill_value=0)
        st_site = stp[stp['siteID'] == site]

        for year in years:
            year_keys = [(s, p, y) for s, p, y in site_keys if y == year]
            # Load spectra
            plot_spectra = {}
            for s, p, y in year_keys:
                spec = load_plot_mean_spectrum(all_clips[(s, p, y)])
                if spec is not None:
                    plot_spectra[p] = spec

            common = sorted(set(plot_spectra.keys()) & set(pivot.index) & set(st_site['plotID']))
            if len(common) < 5:
                continue

            # Distances
            abund = pivot.loc[common].values.astype(float)
            tax_d = pdist(abund, 'braycurtis')

            # Ensure all spectra have same length (trim to minimum)
            spec_list = [plot_spectra[p] for p in common]
            min_bands = min(len(s) for s in spec_list)
            sp = np.array([s[:min_bands] for s in spec_list])
            n_pc = min(5, len(common)-1, sp.shape[1])
            sp_pca = PCA(n_components=n_pc).fit_transform(StandardScaler().fit_transform(sp))
            spec_d = pdist(sp_pca)

            st_data = st_site[st_site['plotID'].isin(common)].set_index('plotID').loc[common][sv].values
            struct_d = pdist(StandardScaler().fit_transform(st_data))

            r_spec, p_spec = pearsonr(tax_d, spec_d)
            r_struct, p_struct = pearsonr(tax_d, struct_d)

            results.append({
                'site': site, 'year': year, 'n_plots': len(common),
                'r_spec_tax': r_spec, 'p_spec_tax': p_spec,
                'r_struct_tax': r_struct, 'p_struct_tax': p_struct,
                'source': 'brdf_001' if year <= 2021 else 'neon_002',
            })

        if results and results[-1]['site'] == site:
            latest = [r for r in results if r['site'] == site]
            print(f"  {site}: {len(latest)} years, "
                  f"spec~tax mean r={np.mean([r['r_spec_tax'] for r in latest]):.3f}, "
                  f"struct~tax mean r={np.mean([r['r_struct_tax'] for r in latest]):.3f}")

    if results:
        rdf = pd.DataFrame(results)
        rdf.to_csv(OUT_DIR / "timeseries_pairwise_beta.csv", index=False)
        print(f"\n  Total: {len(rdf)} site-years")
        print(f"  Spec~Tax:   mean r={rdf['r_spec_tax'].mean():.3f}, sig {(rdf['p_spec_tax']<0.05).sum()}/{len(rdf)}")
        print(f"  Struct~Tax: mean r={rdf['r_struct_tax'].mean():.3f}, sig {(rdf['p_struct_tax']<0.05).sum()}/{len(rdf)}")


# ═══════════════════════════════════════════════════════════════
# 3) Disturbance + spectral change
# ═══════════════════════════════════════════════════════════════

def analysis_3_disturbance_spectral():
    """How does disturbance affect spectral-taxonomic relationship?"""
    print("\n=== Analysis 3: Disturbance x spectral beta ===")

    ts = pd.read_csv(OUT_DIR / "timeseries_pairwise_beta.csv")
    mort = pd.read_csv(OUT_DIR / "plot_mortality_timeseries.csv")

    # Per site-year: get mean mortality
    site_mort = mort.groupby(['siteID', 'year']).agg(
        mean_mort=('mortality_rate', 'mean'),
        max_mort=('mortality_rate', 'max'),
    ).reset_index()

    ts = ts.merge(site_mort.rename(columns={'siteID': 'site'}), on=['site', 'year'], how='left')
    ts['mort_class'] = pd.cut(ts['mean_mort'].fillna(0),
                               bins=[-0.01, 0.05, 0.15, 1.01],
                               labels=['Low', 'Medium', 'High'])

    print(f"  Site-years with mortality data: {ts['mean_mort'].notna().sum()}/{len(ts)}")
    print(ts['mort_class'].value_counts())

    # Compare r_spec_tax by mortality class
    print("\n  Spectral~Taxonomic r by mortality:")
    for cls in ['Low', 'Medium', 'High']:
        sub = ts[ts['mort_class'] == cls]
        if len(sub) > 0:
            print(f"    {cls}: mean r={sub['r_spec_tax'].mean():.3f} (n={len(sub)})")

    print("\n  Structural~Taxonomic r by mortality:")
    for cls in ['Low', 'Medium', 'High']:
        sub = ts[ts['mort_class'] == cls]
        if len(sub) > 0:
            print(f"    {cls}: mean r={sub['r_struct_tax'].mean():.3f} (n={len(sub)})")

    ts.to_csv(OUT_DIR / "disturbance_spectral_beta.csv", index=False)

    # ── Figure ──
    fig, axes = plt.subplots(1, 3, figsize=(21, 7))

    # Panel 1: r_struct vs r_spec across years
    ax = axes[0]
    for site in sorted(ts['site'].unique()):
        sub = ts[ts['site'] == site].sort_values('year')
        ax.plot(sub['year'], sub['r_struct_tax'], 'o-', alpha=0.5, markersize=3)
    ax.axhline(y=0, color='black', linewidth=0.5, linestyle='--')
    ax.set_xlabel('Year')
    ax.set_ylabel('Mantel r')
    ax.set_title('Structural ~ Taxonomic beta\nacross years', fontweight='bold')

    ax2 = axes[1]
    for site in sorted(ts['site'].unique()):
        sub = ts[ts['site'] == site].sort_values('year')
        ax2.plot(sub['year'], sub['r_spec_tax'], 's-', alpha=0.5, markersize=3)
    ax2.axhline(y=0, color='black', linewidth=0.5, linestyle='--')
    ax2.set_xlabel('Year')
    ax2.set_ylabel('Mantel r')
    ax2.set_title('Spectral ~ Taxonomic beta\nacross years (BRDF corrected)', fontweight='bold')

    # Panel 3: by mortality class
    ax3 = axes[2]
    colors = {'Low': '#2ca02c', 'Medium': '#ff7f0e', 'High': '#d62728'}
    for cls in ['Low', 'Medium', 'High']:
        sub = ts[ts['mort_class'] == cls]
        if sub.empty: continue
        ax3.scatter(sub['r_struct_tax'], sub['r_spec_tax'], c=colors[cls],
                   s=30, alpha=0.6, label=f'{cls} mortality (n={len(sub)})')
    ax3.axhline(y=0, color='gray', linewidth=0.5)
    ax3.axvline(x=0, color='gray', linewidth=0.5)
    ax3.set_xlabel('Structural ~ Taxonomic r')
    ax3.set_ylabel('Spectral ~ Taxonomic r')
    ax3.set_title('Structural vs Spectral predictive power\nby disturbance level', fontweight='bold')
    ax3.legend(fontsize=9)

    plt.suptitle('Time-series Pairwise Beta Diversity (BRDF-corrected .001 + NEON .002)\n19 sites, 2013-2025',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(FIG_DIR / 'brdf_timeseries_beta.png', dpi=200, bbox_inches='tight')
    plt.savefig('C:/Users/star1/Documents/GitHub/NEON_Resilience/docs/brdf_timeseries_beta.png',
                dpi=200, bbox_inches='tight')
    print("\n  Figure saved!")


if __name__ == "__main__":
    analysis_1_brdf_comparison()
    analysis_2_timeseries_beta()
    analysis_3_disturbance_spectral()
    print("\nAll analyses complete.")
