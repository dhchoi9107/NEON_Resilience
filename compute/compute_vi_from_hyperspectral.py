"""
Compute Vegetation Indices from BRDF-corrected Hyperspectral (426 bands)
========================================================================
Reads plot-level 100x100m hyperspectral clips, computes NDVI, EVI, ARVI,
PRI, SAVI, LAI(proxy), fPAR(proxy) at 1m resolution, and outputs in the
same format as plot_spectral_1m.csv.

This fills the 2013-2021 gap where NEON VI products weren't downloaded.

Sources:
  - E:/neon_lidar/hyperspectral_brdf_corrected/{SITE}/ (2013-2021)
  - E:/neon_lidar/hyperspectral_plots/{SITE}/         (2022+)

Usage:
  python compute/compute_vi_from_hyperspectral.py
  python compute/compute_vi_from_hyperspectral.py --site HARV
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import argparse
import re
import warnings

import numpy as np
import pandas as pd
import rasterio
from scipy.spatial import ConvexHull
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

from site_config import SITES, SPEC_DIV_DIR

warnings.filterwarnings("ignore")

BRDF_DIR = Path("E:/neon_lidar/hyperspectral_brdf_corrected")
PLOTS_DIR = Path("E:/neon_lidar/hyperspectral_plots")
OUT_PATH = SPEC_DIV_DIR / "plot_spectral_1m.csv"

PLOT_HALF = 20  # 40m plot = ±20m from center
GRAIN_SIZES = [1, 2, 5, 10]
MAX_SAMPLE_RAO = 500
MAX_SAMPLE_PCA = 2000

# NEON AOP: 426 bands, ~380-2510nm, ~5nm spacing
WL_START = 380
WL_STEP = 5
N_BANDS = 426

WAVELENGTHS = np.arange(WL_START, WL_START + N_BANDS * WL_STEP, WL_STEP)[:N_BANDS]


def find_band_idx(target_nm):
    """Find 0-indexed band number for target wavelength."""
    return int(np.argmin(np.abs(WAVELENGTHS - target_nm)))


# Pre-compute band indices for VI calculations
B_BLUE = find_band_idx(475)
B_RED = find_band_idx(670)
B_NIR = find_band_idx(850)
B_GREEN = find_band_idx(555)
B_531 = find_band_idx(531)
B_570 = find_band_idx(570)
B_RE = find_band_idx(720)   # red-edge for LAI proxy


def compute_vis(data_cube):
    """Compute VIs from (bands, H, W) hyperspectral cube.

    Returns dict of 2D arrays: {vi_name: (H, W) array}
    Data is int16 reflectance * 10000 scale.
    """
    # Convert to float reflectance [0, 1]
    cube = data_cube.astype(np.float32) / 10000.0

    blue = cube[B_BLUE]
    red = cube[B_RED]
    nir = cube[B_NIR]
    green = cube[B_GREEN]
    r531 = cube[B_531]
    r570 = cube[B_570]
    r_re = cube[B_RE]

    # Mask invalid pixels
    valid = (data_cube[B_RED] > 0) & (data_cube[B_RED] < 10000) & \
            (data_cube[B_NIR] > 0) & (data_cube[B_NIR] < 10000)

    def safe_ratio(a, b):
        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.where(np.abs(b) > 1e-8, a / b, np.nan)
        result[~valid] = np.nan
        return result

    vis = {}

    # NDVI = (NIR - Red) / (NIR + Red)
    vis["NDVI"] = safe_ratio(nir - red, nir + red)

    # EVI = 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)
    evi_denom = nir + 6.0 * red - 7.5 * blue + 1.0
    vis["EVI"] = np.where(valid, 2.5 * safe_ratio(nir - red, evi_denom), np.nan)

    # ARVI = (NIR - (2*Red - Blue)) / (NIR + (2*Red - Blue))
    rb = 2.0 * red - blue
    vis["ARVI"] = safe_ratio(nir - rb, nir + rb)

    # PRI = (R531 - R570) / (R531 + R570)
    vis["PRI"] = safe_ratio(r531 - r570, r531 + r570)

    # SAVI = 1.5 * (NIR - Red) / (NIR + Red + 0.5)
    vis["SAVI"] = np.where(valid, 1.5 * safe_ratio(nir - red, nir + red + 0.5), np.nan)

    # LAI proxy: using NDVI-based empirical relationship
    # LAI ~ -1/k * ln((NDVI_inf - NDVI) / (NDVI_inf - NDVI_soil))
    # Simplified: LAI ~ 3.618 * EVI - 0.118 (Boegh et al. 2002)
    vis["LAI"] = np.where(valid, 3.618 * vis["EVI"] - 0.118, np.nan)

    # fPAR proxy: fPAR ~ 1.24 * NDVI - 0.168 (Myneni & Williams 1994)
    vis["fPAR"] = np.where(valid,
                           np.clip(1.24 * vis["NDVI"] - 0.168, 0, 1),
                           np.nan)

    return vis


# ── Diversity metrics (same as compute_plot_spectral_1m.py) ──────────────

def compute_rao_q(data, n_sample=MAX_SAMPLE_RAO):
    n = len(data)
    if n < 2:
        return np.nan
    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        data = data[idx]
    diff = data[:, np.newaxis, :] - data[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=2))
    mask = np.triu(np.ones(dists.shape, dtype=bool), k=1)
    return float(np.mean(dists[mask]))


def compute_spectral_cv(data):
    if len(data) < 2:
        return np.nan
    means = np.mean(data, axis=0)
    sds = np.std(data, axis=0)
    cvs = np.where(np.abs(means) > 1e-6, sds / np.abs(means), np.nan)
    return float(np.nanmean(cvs))


def compute_spectral_shannon(data, n_bins=20):
    if len(data) < 2:
        return np.nan
    ents = []
    for b in range(data.shape[1]):
        vals = data[:, b]
        counts, _ = np.histogram(vals, bins=n_bins)
        p = counts / counts.sum()
        p = p[p > 0]
        ents.append(-np.sum(p * np.log(p)))
    return float(np.mean(ents))


def compute_pca_diversity(data, n_pcs=3):
    result = {"spectral_FRic": np.nan, "spectral_FDiv": np.nan, "spectral_FEve": np.nan}
    n = len(data)
    if n < n_pcs + 2:
        return result

    sample = data
    if n > MAX_SAMPLE_PCA:
        idx = np.random.choice(n, MAX_SAMPLE_PCA, replace=False)
        sample = data[idx]

    pca = PCA(n_components=min(n_pcs, sample.shape[1]))
    scores = pca.fit_transform(sample)

    try:
        if scores.shape[1] >= 3 and len(scores) >= 4:
            hull = ConvexHull(scores[:, :3])
            result["spectral_FRic"] = float(hull.volume)
        elif scores.shape[1] >= 2 and len(scores) >= 3:
            hull = ConvexHull(scores[:, :2])
            result["spectral_FRic"] = float(hull.volume)
    except Exception:
        pass

    centroid = np.mean(scores, axis=0)
    dists = np.sqrt(np.sum((scores - centroid) ** 2, axis=1))
    if np.max(dists) > 0:
        result["spectral_FDiv"] = float(np.mean(dists) / np.max(dists))

    if len(scores) >= 3:
        sc = scores[:min(len(scores), MAX_SAMPLE_PCA)]
        dm = squareform(pdist(sc))
        mst = minimum_spanning_tree(dm)
        edges = mst.data[mst.data > 0]
        if len(edges) > 1:
            total = edges.sum()
            pev = edges / total
            expected = 1.0 / len(edges)
            eve = 1.0 - np.sum(np.abs(pev - expected)) / (2.0 * (1.0 - expected))
            result["spectral_FEve"] = float(np.clip(eve, 0, 1))

    return result


# ── Processing ───────────────────────────────────────────────────────────

VI_NAMES = ["NDVI", "EVI", "ARVI", "PRI", "SAVI", "LAI", "fPAR"]


def process_plot_file(tif_path, site, plot_id, year, grain_sizes):
    """Process one hyperspectral plot clip → VI metrics at each grain."""
    with rasterio.open(str(tif_path)) as src:
        n_bands = src.count
        h, w = src.height, src.width
        data_cube = src.read()  # (426, H, W)

    if n_bands != N_BANDS:
        return []

    # Extract 40x40m center
    half = min(PLOT_HALF, h // 2, w // 2)
    if half < 5:
        return []
    cx, cy = w // 2, h // 2
    center = data_cube[:, cy - half:cy + half, cx - half:cx + half]

    # Compute VIs
    vis = compute_vis(center)

    # Stack VIs into (n_vi, H, W)
    vi_stack = np.stack([vis[name] for name in VI_NAMES], axis=0)

    results = []
    for grain in grain_sizes:
        # Aggregate to grain
        if grain == 1:
            n_vi, gh, gw = vi_stack.shape
            flat = vi_stack.reshape(n_vi, -1).T  # (n_pixels, n_vi)
            valid = ~np.any(np.isnan(flat), axis=1)
            flat = flat[valid]
        else:
            n_vi, gh, gw = vi_stack.shape
            bh = (gh // grain) * grain
            bw = (gw // grain) * grain
            trimmed = vi_stack[:, :bh, :bw]
            reshaped = trimmed.reshape(n_vi, bh // grain, grain, bw // grain, grain)
            with np.errstate(all='ignore'):
                agg = np.nanmean(reshaped, axis=(2, 4))
            flat = agg.reshape(n_vi, -1).T
            valid = ~np.any(np.isnan(flat), axis=1)
            flat = flat[valid]

        if len(flat) < 4:
            continue

        row = {
            "siteID": site,
            "plotID": plot_id,
            "year": year,
            "grain_m": grain,
            "n_pixels": len(flat),
        }

        # Per-band stats
        for i, name in enumerate(VI_NAMES):
            row[f"{name}_mean"] = float(np.nanmean(flat[:, i]))
            row[f"{name}_sd"] = float(np.nanstd(flat[:, i]))

        # Diversity metrics
        row["rao_q"] = compute_rao_q(flat)
        row["spectral_cv"] = compute_spectral_cv(flat)
        row["spectral_shannon"] = compute_spectral_shannon(flat)

        pca_m = compute_pca_diversity(flat)
        row.update(pca_m)

        results.append(row)

    return results


def process_site(site, grain_sizes):
    """Process all hyperspectral clips for a site from both directories."""
    all_rows = []

    # Source 1: BRDF corrected (2013-2021)
    brdf_dir = BRDF_DIR / site
    if brdf_dir.exists():
        tifs = sorted(brdf_dir.glob("*_hyper_brdf.tif"))
        by_year = {}
        for f in tifs:
            m = re.match(r"(\d{4})_(.+)_hyper_brdf\.tif", f.name)
            if m:
                yr, pid = int(m.group(1)), m.group(2)
                by_year.setdefault(yr, []).append((pid, f))

        for yr in sorted(by_year):
            n_ok = 0
            for pid, fpath in by_year[yr]:
                rows = process_plot_file(fpath, site, pid, yr, grain_sizes)
                all_rows.extend(rows)
                if rows:
                    n_ok += 1
            if n_ok > 0:
                print(f"    BRDF {yr}: {n_ok}/{len(by_year[yr])} plots")

    # Source 2: hyperspectral_plots (2022+, non-BRDF clips)
    plots_dir = PLOTS_DIR / site
    if plots_dir.exists():
        tifs = sorted(plots_dir.glob("*_hyper.tif"))
        by_year = {}
        for f in tifs:
            m = re.match(r"(\d{4})_(.+)_hyper\.tif", f.name)
            if m:
                yr, pid = int(m.group(1)), m.group(2)
                by_year.setdefault(yr, []).append((pid, f))

        for yr in sorted(by_year):
            n_ok = 0
            for pid, fpath in by_year[yr]:
                rows = process_plot_file(fpath, site, pid, yr, grain_sizes)
                all_rows.extend(rows)
                if rows:
                    n_ok += 1
            if n_ok > 0:
                print(f"    Plots {yr}: {n_ok}/{len(by_year[yr])} plots")

    return all_rows


def main():
    parser = argparse.ArgumentParser(
        description="Compute VIs from hyperspectral for plot_spectral_1m")
    parser.add_argument("--site", type=str, default=None)
    parser.add_argument("--grain", type=int, default=None)
    args = parser.parse_args()

    SPEC_DIV_DIR.mkdir(parents=True, exist_ok=True)

    sites = [args.site] if args.site else SITES
    grains = [args.grain] if args.grain else GRAIN_SIZES

    print(f"Computing VIs from Hyperspectral Data")
    print(f"  BRDF dir: {BRDF_DIR}")
    print(f"  Plots dir: {PLOTS_DIR}")
    print(f"  Grains: {grains}")
    print()

    all_results = []
    for site in sites:
        print(f"  {site}:")
        rows = process_site(site, grains)
        all_results.extend(rows)
        n1m = sum(1 for r in rows if r["grain_m"] == 1)
        if n1m:
            print(f"    → {n1m} plot-years (1m)")

    if not all_results:
        print("No results!")
        return

    hyper_df = pd.DataFrame(all_results)

    # Load existing VI-based spectral data
    existing_path = OUT_PATH
    if existing_path.exists():
        existing = pd.read_csv(existing_path)
        print(f"\nExisting plot_spectral_1m.csv: {len(existing)} rows")

        # Remove any overlapping site-year-plot-grain combos
        hyper_keys = set(zip(hyper_df["siteID"], hyper_df["plotID"],
                             hyper_df["year"], hyper_df["grain_m"]))
        existing_mask = ~existing.apply(
            lambda r: (r["siteID"], r["plotID"], int(r["year"]), int(r["grain_m"])) in hyper_keys,
            axis=1
        )
        existing_kept = existing[existing_mask]
        print(f"  Kept (no overlap): {len(existing_kept)} rows")
        print(f"  New from hyperspectral: {len(hyper_df)} rows")

        # Align columns - add missing columns as NaN
        for col in existing.columns:
            if col not in hyper_df.columns:
                hyper_df[col] = np.nan
        for col in hyper_df.columns:
            if col not in existing_kept.columns:
                existing_kept[col] = np.nan

        combined = pd.concat([existing_kept, hyper_df], ignore_index=True)
    else:
        combined = hyper_df

    # Sort
    combined = combined.sort_values(["siteID", "plotID", "year", "grain_m"]).reset_index(drop=True)

    # Backup (replace existing backup on Windows)
    if existing_path.exists():
        backup = existing_path.with_suffix(".csv.bak")
        if backup.exists():
            backup.unlink()
        existing_path.rename(backup)
        print(f"  Backed up to {backup.name}")

    combined.to_csv(OUT_PATH, index=False)
    print(f"\nSaved: {OUT_PATH} ({len(combined)} rows)")

    # Summary
    g1 = combined[combined["grain_m"] == 1]
    print(f"\n1m grain summary:")
    print(f"  Total plot-years: {len(g1)}")
    print(f"  Sites: {g1['siteID'].nunique()}")
    print(f"  Years: {sorted(g1['year'].unique())}")
    print(f"  Site-years: {g1.groupby(['siteID','year']).ngroups}")


if __name__ == "__main__":
    main()
