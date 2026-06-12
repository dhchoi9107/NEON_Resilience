"""
Hyperspectral Spectral Diversity (426 bands, 1m, plot-level)
==============================================================
Computes spectral diversity from NEON BRDF-corrected hyperspectral
reflectance (DP3.30006.002) clipped to plot footprints (100x100m).

Processing steps:
  1. Bad band removal (water absorption: 1350-1460nm, 1790-2000nm, >2400nm)
  2. Nodata masking (reflectance <= 0 or >= 10000)
  3. PCA (site-level or pooled)
  4. Per-plot diversity metrics at multiple grain sizes

Metrics:
  - PCA-based: FRic, FDiv, FEve (Schneider et al. 2017)
  - Rao's Q on PC space (Rocchini et al. 2017)
  - Spectral alpha/beta/gamma partitioning (Laliberte et al. 2020)
  - Functional trait indices (chlorophyll, water, nitrogen proxies)

Usage:
  python compute/compute_hyperspectral_diversity.py
  python compute/compute_hyperspectral_diversity.py --site HARV
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import warnings
import re

import numpy as np
import pandas as pd
import rasterio
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform

from site_config import SITES, SPEC_DIV_DIR

warnings.filterwarnings("ignore")

HYPER_DIR = Path("E:/neon_lidar/hyperspectral_plots")
PLOT_HALF = 20  # 40m plot half-width (center of 100x100m clip)
GRAIN_SIZES = [1, 2, 5, 10]
MAX_SAMPLE = 2000
N_PCS = 5  # keep more PCs for richer trait space

# NEON AOP wavelengths: 426 bands from ~380nm to ~2510nm, ~5nm spacing
# Bad bands: water absorption + noisy edges
BAD_RANGES_NM = [
    (0, 395),       # noisy UV edge
    (1340, 1465),   # water absorption
    (1790, 2005),   # water absorption
    (2400, 2600),   # noisy SWIR edge
]


# ─── Preprocessing ──────────────────────────────────────────────────────────

def get_good_band_mask(n_bands=426, wl_start=380, wl_step=5):
    """Return boolean mask for good (non-water-absorption) bands."""
    wavelengths = np.arange(wl_start, wl_start + n_bands * wl_step, wl_step)[:n_bands]
    good = np.ones(n_bands, dtype=bool)
    for lo, hi in BAD_RANGES_NM:
        good &= ~((wavelengths >= lo) & (wavelengths <= hi))
    return good, wavelengths


def load_plot_hyper(tif_path, grain_m=1):
    """Load hyperspectral plot clip, remove bad bands, mask nodata.

    Returns (n_pixels, n_good_bands) array or None.
    """
    with rasterio.open(str(tif_path)) as src:
        n_bands = src.count
        # Read center 40x40m (plot footprint) from clip
        # Clip may be smaller than 100x100 if plot is at tile edge
        w, h = src.width, src.height
        half = min(PLOT_HALF, w // 2, h // 2)
        if half < 5:
            return None, None
        cx, cy = w // 2, h // 2
        window = rasterio.windows.Window(
            cx - half, cy - half, 2 * half, 2 * half
        )
        data = src.read(window=window).astype(np.float32)  # (bands, H, W)

    good_mask, wavelengths = get_good_band_mask(n_bands)
    data = data[good_mask]  # (n_good_bands, H, W)
    n_good = data.shape[0]

    # Nodata mask: reflectance should be 0-10000 (NEON scale factor)
    # Values <= 0 or >= 10000 are nodata/shadow/saturated
    valid_pixel = np.all((data > 0) & (data < 10000), axis=0)  # (H, W)

    # Aggregate to grain size
    if grain_m > 1:
        h, w = data.shape[1], data.shape[2]
        bh = (h // grain_m) * grain_m
        bw = (w // grain_m) * grain_m
        data = data[:, :bh, :bw]
        valid_pixel = valid_pixel[:bh, :bw]

        # Reshape to blocks
        data = data.reshape(n_good, bh // grain_m, grain_m, bw // grain_m, grain_m)
        data = np.nanmean(data, axis=(2, 4))

        valid_pixel = valid_pixel.reshape(bh // grain_m, grain_m, bw // grain_m, grain_m)
        valid_pixel = np.all(valid_pixel, axis=(1, 3))

    # Flatten to (n_pixels, n_bands)
    flat = data.reshape(n_good, -1).T  # (n_pixels, n_bands)
    valid_flat = valid_pixel.ravel()
    flat = flat[valid_flat]

    if len(flat) < 10:
        return None, wavelengths[good_mask] if len(wavelengths) > 0 else None

    return flat, wavelengths[good_mask]


# ─── Diversity metrics ───────────────────────────────────────────────────────

def compute_pca_diversity(data, n_pcs=N_PCS):
    """PCA-based functional spectral diversity."""
    result = {"hyper_FRic": np.nan, "hyper_FDiv": np.nan, "hyper_FEve": np.nan,
              "hyper_var_explained": np.nan}

    n = len(data)
    if n < n_pcs + 2:
        return result, None

    sample = data
    if n > MAX_SAMPLE:
        idx = np.random.choice(n, MAX_SAMPLE, replace=False)
        sample = data[idx]

    # Standardize before PCA (reflectance scale 0-10000 → z-score)
    from sklearn.preprocessing import StandardScaler
    sample = StandardScaler().fit_transform(sample)

    pca = PCA(n_components=min(n_pcs, sample.shape[1]))
    scores = pca.fit_transform(sample)
    result["hyper_var_explained"] = float(np.sum(pca.explained_variance_ratio_))

    # FRic: convex hull volume in first 3 PCs
    if scores.shape[1] >= 3 and len(scores) >= 4:
        try:
            hull = ConvexHull(scores[:, :3])
            result["hyper_FRic"] = float(hull.volume)
        except Exception:
            pass

    # FDiv
    centroid = np.mean(scores, axis=0)
    dists = np.sqrt(np.sum((scores - centroid) ** 2, axis=1))
    if np.max(dists) > 0:
        result["hyper_FDiv"] = float(np.mean(dists) / np.max(dists))

    # FEve
    if len(scores) >= 3:
        sc = scores[:min(len(scores), MAX_SAMPLE)]
        dm = squareform(pdist(sc))
        mst = minimum_spanning_tree(dm)
        edges = mst.data[mst.data > 0]
        if len(edges) > 1:
            total = edges.sum()
            pev = edges / total
            expected = 1.0 / len(edges)
            eve = 1.0 - np.sum(np.abs(pev - expected)) / (2.0 * (1.0 - expected))
            result["hyper_FEve"] = float(np.clip(eve, 0, 1))

    return result, scores


def compute_rao_q(data, n_sample=500):
    n = len(data)
    if n < 2:
        return np.nan
    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        data = data[idx]
    dists = pdist(data, metric="euclidean")
    return float(np.mean(dists))


def compute_spectral_cv(data):
    if len(data) < 2:
        return np.nan
    means = np.mean(data, axis=0)
    sds = np.std(data, axis=0)
    cvs = np.where(np.abs(means) > 1e-6, sds / np.abs(means), np.nan)
    return float(np.nanmean(cvs))


def compute_trait_indices(data, wavelengths):
    """Compute functional trait proxy indices from hyperspectral data.

    Returns dict with per-plot mean and SD of each index.
    """
    result = {}
    if wavelengths is None or len(data) == 0:
        return result

    wl = wavelengths

    def band_idx(target_nm):
        return np.argmin(np.abs(wl - target_nm))

    def safe_ratio(a, b):
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(np.abs(b) > 1e-6, a / b, np.nan)
        return r

    try:
        # NDVI: (R850 - R670) / (R850 + R670)
        r670 = data[:, band_idx(670)]
        r850 = data[:, band_idx(850)]
        ndvi = safe_ratio(r850 - r670, r850 + r670)
        result["trait_NDVI_mean"] = float(np.nanmean(ndvi))
        result["trait_NDVI_sd"] = float(np.nanstd(ndvi))

        # NDWI: (R860 - R1240) / (R860 + R1240)
        r860 = data[:, band_idx(860)]
        r1240 = data[:, band_idx(1240)]
        ndwi = safe_ratio(r860 - r1240, r860 + r1240)
        result["trait_NDWI_mean"] = float(np.nanmean(ndwi))
        result["trait_NDWI_sd"] = float(np.nanstd(ndwi))

        # CHL index: R750/R710 (chlorophyll)
        r750 = data[:, band_idx(750)]
        r710 = data[:, band_idx(710)]
        chl = safe_ratio(r750, r710)
        result["trait_CHL_mean"] = float(np.nanmean(chl))
        result["trait_CHL_sd"] = float(np.nanstd(chl))

        # NDNI: (R1510 - R1680) / (R1510 + R1680) — nitrogen
        r1510 = data[:, band_idx(1510)]
        r1680 = data[:, band_idx(1680)]
        ndni = safe_ratio(np.log(1/r1510) - np.log(1/r1680),
                          np.log(1/r1510) + np.log(1/r1680))
        result["trait_NDNI_mean"] = float(np.nanmean(ndni))
        result["trait_NDNI_sd"] = float(np.nanstd(ndni))

        # PRI: (R531 - R570) / (R531 + R570) — photosynthetic efficiency
        r531 = data[:, band_idx(531)]
        r570 = data[:, band_idx(570)]
        pri = safe_ratio(r531 - r570, r531 + r570)
        result["trait_PRI_mean"] = float(np.nanmean(pri))
        result["trait_PRI_sd"] = float(np.nanstd(pri))

    except Exception:
        pass

    return result


# ─── Laliberte partitioning ──────────────────────────────────────────────────

def laliberte_partition(plot_scores_list):
    """Partition spectral diversity into alpha, beta, gamma (Laliberte et al. 2020).

    plot_scores_list: list of (n_pixels, n_pcs) arrays per plot.
    Returns dict with alpha_mean, beta, gamma.
    """
    result = {"spec_alpha": np.nan, "spec_beta": np.nan, "spec_gamma": np.nan}

    valid = [s for s in plot_scores_list if s is not None and len(s) >= 10]
    if len(valid) < 2:
        return result

    # Alpha: mean within-plot Rao's Q
    alphas = [compute_rao_q(s, n_sample=300) for s in valid]
    alphas = [a for a in alphas if not np.isnan(a)]
    if not alphas:
        return result
    result["spec_alpha"] = float(np.mean(alphas))

    # Gamma: Rao's Q of all plots pooled
    pooled = np.vstack([s[:min(len(s), 200)] for s in valid])
    result["spec_gamma"] = compute_rao_q(pooled, n_sample=1000)

    # Beta: gamma - alpha
    result["spec_beta"] = result["spec_gamma"] - result["spec_alpha"]

    return result


# ─── Main processing ────────────────────────────────────────────────────────

def process_site(site, grain_sizes=GRAIN_SIZES):
    """Process all hyperspectral plot clips for a site."""
    site_dir = HYPER_DIR / site
    if not site_dir.exists():
        print(f"  SKIP {site}: no hyperspectral data")
        return []

    tif_files = sorted(site_dir.glob("*_hyper.tif"))
    if not tif_files:
        print(f"  SKIP {site}: no clipped TIFs")
        return []

    # Group by year
    by_year = {}
    for f in tif_files:
        m = re.match(r"(\d{4})_(.+)_hyper\.tif", f.name)
        if not m:
            continue
        yr, plot_id = m.group(1), m.group(2)
        if yr not in by_year:
            by_year[yr] = []
        by_year[yr].append((plot_id, f))

    print(f"  {site}: {len(tif_files)} files, {len(by_year)} years")

    results = []
    for yr in sorted(by_year.keys()):
        plots = by_year[yr]

        # Collect PCA scores for Laliberte partitioning (1m grain)
        plot_scores_1m = []

        for grain in grain_sizes:
            for plot_id, tif_path in plots:
                data, wavelengths = load_plot_hyper(tif_path, grain_m=grain)
                if data is None:
                    continue

                row = {
                    "siteID": site,
                    "plotID": plot_id,
                    "year": int(yr),
                    "grain_m": grain,
                    "n_pixels": len(data),
                    "n_bands": data.shape[1],
                }

                # PCA diversity
                pca_metrics, scores = compute_pca_diversity(data)
                row.update(pca_metrics)

                # Rao's Q on PCA scores (standardized space)
                row["hyper_RaoQ"] = compute_rao_q(scores, n_sample=500) if scores is not None else np.nan

                # Spectral CV
                row["hyper_CV"] = compute_spectral_cv(data)

                # Trait indices (only at 1m — they're per-pixel)
                if grain == 1:
                    traits = compute_trait_indices(data, wavelengths)
                    row.update(traits)

                    # Save PCA scores for Laliberte partitioning
                    if scores is not None:
                        plot_scores_1m.append(scores)

                results.append(row)

        # Laliberte alpha/beta/gamma partitioning (1m only)
        if plot_scores_1m:
            partition = laliberte_partition(plot_scores_1m)
            # Append to all 1m rows for this year
            for r in results:
                if r["siteID"] == site and r["year"] == int(yr) and r["grain_m"] == 1:
                    r.update(partition)

        n_plots = len(plots)
        n_ok = sum(1 for r in results
                   if r["siteID"] == site and r["year"] == int(yr) and r["grain_m"] == 1)
        if n_ok > 0:
            mean_fric = np.nanmean([r["hyper_FRic"] for r in results
                                    if r["siteID"] == site and r["year"] == int(yr)
                                    and r["grain_m"] == 1])
            print(f"    {yr}: {n_ok}/{n_plots} plots, FRic={mean_fric:.1f}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Compute spectral diversity from 426-band hyperspectral data")
    parser.add_argument("--site", type=str, default=None)
    parser.add_argument("--grain", type=int, default=None)
    args = parser.parse_args()

    SPEC_DIV_DIR.mkdir(parents=True, exist_ok=True)

    sites = [args.site] if args.site else SITES
    grains = [args.grain] if args.grain else GRAIN_SIZES

    print(f"Hyperspectral Diversity (426 bands, BRDF corrected)")
    print(f"Input: {HYPER_DIR}")
    print(f"Grains: {grains}")
    print(f"PCs: {N_PCS}\n")

    all_results = []
    for site in sites:
        print(f"\n{'='*50}")
        rows = process_site(site, grain_sizes=grains)
        all_results.extend(rows)

    if all_results:
        df = pd.DataFrame(all_results)
        out_path = SPEC_DIV_DIR / "hyperspectral_diversity.csv"
        df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path} ({len(df)} rows)")

        # Summary
        for g in grains:
            sub = df[df["grain_m"] == g]
            if not sub.empty:
                print(f"  {g}m: {len(sub)} plot-years, "
                      f"FRic={sub['hyper_FRic'].mean():.1f}, "
                      f"RaoQ={sub['hyper_RaoQ'].mean():.0f}, "
                      f"n_px={sub['n_pixels'].mean():.0f}")


if __name__ == "__main__":
    main()
