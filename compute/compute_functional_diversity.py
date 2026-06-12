"""
Functional Diversity: Spectral + Structural Trait Space (1m resolution)
========================================================================
Combines 1m structural diversity (from plot-level FSD rasters) and
1m spectral diversity (from VI products) into a unified functional
trait space per NEON plot.

Uses plot-level 1m FSD rasters (140x140m, from compute_plot_fsd_1m.py)
clipped to the 40x40m plot footprint = ~1600 pixels at 1m.

Traits used:
  Structural (from 1m FSD): mean_max_canopy_ht, FHD, LAI, vert_sd, GC
  Spectral (from plot_spectral_1m.csv): per-plot PCA scores or VI means

Metrics per plot:
  FRic  - Functional richness (convex hull volume)
  FDiv  - Functional divergence
  FEve  - Functional evenness (MST regularity)
  Rao_Q - Mean pairwise trait distance

Usage:
  python compute/compute_functional_diversity.py
  python compute/compute_functional_diversity.py --site HARV
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import warnings

import numpy as np
import pandas as pd
import rasterio
from scipy.spatial import ConvexHull
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform

from site_config import (
    SITES, FUNC_DIV_DIR, VEG_STRUCT_DIR, FSD_BANDS, get_fsd_band_index,
)

warnings.filterwarnings("ignore")

FSD_1M_DIR = Path("E:/neon_lidar/structural_diversity_1m_plots")
PLOT_HALF = 20  # 40m plot half-width
WINDOW_HALF = 70  # FSD 1m raster is plot center ± 70m

# Structural traits to extract from 1m FSD
STRUCT_TRAITS = ["mean_max_canopy_ht", "FHD", "LAI", "vert_sd", "GC"]
STRUCT_BAND_IDX = {t: get_fsd_band_index(t) for t in STRUCT_TRAITS}

MAX_SAMPLE = 2000  # subsample for expensive metrics


# ─── Load plot coordinates ──────────────────────────────────────────────────

def load_plot_coords():
    ppy_path = VEG_STRUCT_DIR / "vst_perplotperyear.csv"
    if not ppy_path.exists():
        return pd.DataFrame()
    df = pd.read_csv(ppy_path, low_memory=False)
    if "easting" not in df.columns:
        return pd.DataFrame()
    coords = df[["siteID", "plotID", "easting", "northing"]].dropna()
    return coords.drop_duplicates(subset=["siteID", "plotID"])


# ─── Extract structural traits from 1m FSD plot raster ──────────────────────

def extract_struct_traits_1m(site, year, plot_id, plot_x, plot_y):
    """Extract structural traits from 1m FSD raster within 40x40m plot.

    Returns (n_pixels, n_traits) array or None.
    """
    fsd_path = FSD_1M_DIR / site / f"{year}_{plot_id}_FSD_1m.tif"
    if not fsd_path.exists():
        return None

    with rasterio.open(str(fsd_path)) as src:
        # The raster is 140x140m centered on the plot
        # Extract only the 40x40m plot area (center ±20m)
        xmin = plot_x - PLOT_HALF
        xmax = plot_x + PLOT_HALF
        ymin = plot_y - PLOT_HALF
        ymax = plot_y + PLOT_HALF

        try:
            window = rasterio.windows.from_bounds(
                xmin, ymin, xmax, ymax, transform=src.transform
            )
            window = window.intersection(
                rasterio.windows.Window(0, 0, src.width, src.height)
            )
        except Exception:
            return None

        if window.width < 5 or window.height < 5:
            return None

        columns = []
        for trait in STRUCT_TRAITS:
            bidx = STRUCT_BAND_IDX[trait]
            data = src.read(bidx, window=window).astype(np.float32).ravel()
            columns.append(data)

    matrix = np.column_stack(columns)
    # Remove rows with any NaN
    valid = ~np.any(np.isnan(matrix), axis=1)
    matrix = matrix[valid]

    return matrix if len(matrix) >= 10 else None


# ─── Functional diversity metrics ────────────────────────────────────────────

def compute_fric(matrix):
    n, d = matrix.shape
    if n < d + 1:
        return np.nan
    try:
        dims = min(d, 3)
        hull = ConvexHull(matrix[:, :dims])
        return float(hull.volume)
    except Exception:
        return np.nan


def compute_fdiv(matrix):
    centroid = np.mean(matrix, axis=0)
    dists = np.sqrt(np.sum((matrix - centroid) ** 2, axis=1))
    max_d = np.max(dists)
    if max_d < 1e-8:
        return np.nan
    return float(np.mean(dists) / max_d)


def compute_feve(matrix):
    if len(matrix) < 3:
        return np.nan
    m = matrix
    if len(m) > MAX_SAMPLE:
        idx = np.random.choice(len(m), MAX_SAMPLE, replace=False)
        m = m[idx]
    dm = squareform(pdist(m))
    mst = minimum_spanning_tree(dm)
    edges = mst.data[mst.data > 0]
    if len(edges) < 2:
        return np.nan
    total = edges.sum()
    pev = edges / total
    expected = 1.0 / len(edges)
    eve = 1.0 - np.sum(np.abs(pev - expected)) / (2.0 * (1.0 - expected))
    return float(np.clip(eve, 0, 1))


def compute_rao_q(matrix, n_sample=500):
    n = len(matrix)
    if n < 2:
        return np.nan
    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        matrix = matrix[idx]
    dists = pdist(matrix, metric="euclidean")
    return float(np.mean(dists))


def functional_metrics(matrix):
    # Z-score normalize
    means = np.mean(matrix, axis=0)
    sds = np.std(matrix, axis=0)
    sds[sds < 1e-8] = 1.0
    normed = (matrix - means) / sds

    return {
        "func_FRic": compute_fric(normed),
        "func_FDiv": compute_fdiv(normed),
        "func_FEve": compute_feve(normed),
        "func_RaoQ": compute_rao_q(normed),
        "func_n_pixels": len(matrix),
        "func_n_traits": matrix.shape[1],
    }


# ─── Beta functional diversity ──────────────────────────────────────────────

def functional_beta(plot_matrices):
    centroids = []
    for m in plot_matrices:
        if m is not None and len(m) > 0:
            # Normalize per-plot then take centroid
            means = np.mean(m, axis=0)
            sds = np.std(m, axis=0)
            sds[sds < 1e-8] = 1.0
            normed = (m - means) / sds
            centroids.append(np.mean(normed, axis=0))

    if len(centroids) < 2:
        return {"func_beta_rao": np.nan}

    centroids = np.array(centroids)
    dists = pdist(centroids, metric="euclidean")
    return {"func_beta_rao": float(np.mean(dists))}


# ─── Main processing ────────────────────────────────────────────────────────

def process_site(site, plot_coords):
    site_plots = plot_coords[plot_coords["siteID"] == site]
    if site_plots.empty:
        print(f"  SKIP {site}: no plot coordinates")
        return []

    site_fsd_dir = FSD_1M_DIR / site
    if not site_fsd_dir.exists():
        print(f"  SKIP {site}: no 1m FSD directory")
        return []

    # Find available years from FSD filenames
    fsd_files = list(site_fsd_dir.glob("*_FSD_1m.tif"))
    years = sorted(set(int(f.name.split("_")[0]) for f in fsd_files))
    if not years:
        print(f"  SKIP {site}: no 1m FSD files")
        return []

    print(f"  {site}: {len(site_plots)} plots, {len(years)} years")

    results = []
    for year in years:
        plot_matrices = []
        year_results = []

        for _, prow in site_plots.iterrows():
            px, py = prow["easting"], prow["northing"]
            plot_id = prow["plotID"]

            matrix = extract_struct_traits_1m(site, year, plot_id, px, py)
            if matrix is None:
                continue

            metrics = functional_metrics(matrix)
            metrics["siteID"] = site
            metrics["plotID"] = plot_id
            metrics["year"] = year
            year_results.append(metrics)
            plot_matrices.append(matrix)

        # Site-level beta
        beta = functional_beta(plot_matrices)
        for r in year_results:
            r.update(beta)

        results.extend(year_results)

        if year_results:
            mean_fric = np.nanmean([r["func_FRic"] for r in year_results])
            print(f"    {year}: {len(year_results)}/{len(site_plots)} plots, "
                  f"mean FRic={mean_fric:.2f}, n_px={year_results[0]['func_n_pixels']}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Compute functional diversity from 1m data")
    parser.add_argument("--site", type=str, default=None)
    args = parser.parse_args()

    FUNC_DIV_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading plot coordinates...")
    plot_coords = load_plot_coords()
    if plot_coords.empty:
        print("ERROR: No plot coordinates.")
        return
    print(f"  {len(plot_coords)} plots\n")

    sites = [args.site] if args.site else SITES
    all_results = []

    for site in sites:
        print(f"\n{'='*50}")
        rows = process_site(site, plot_coords)
        all_results.extend(rows)

    if all_results:
        df = pd.DataFrame(all_results)
        csv_path = FUNC_DIV_DIR / "functional_diversity_all.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path} ({len(df)} plot-years)")

        # Summary
        summary = df.groupby("siteID").agg(
            n_plotyears=("plotID", "count"),
            FRic_mean=("func_FRic", "mean"),
            FDiv_mean=("func_FDiv", "mean"),
            RaoQ_mean=("func_RaoQ", "mean"),
            n_pixels_mean=("func_n_pixels", "mean"),
        ).reset_index()
        print("\nSite summary:")
        print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
