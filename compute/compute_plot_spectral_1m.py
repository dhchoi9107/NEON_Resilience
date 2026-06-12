"""
Plot-level Spectral Diversity at 1m Resolution
================================================
Extracts 1m VI pixels within each NEON 40x40m plot footprint and computes
spectral diversity metrics at multiple grain sizes (1m, 2m, 5m, 10m).

This addresses the scale-dependency issue: 10m resolution gives only 16
pixels per plot (insufficient for ConvexHull), while 1m gives 1,600 pixels.

Metrics per plot, per grain size:
  - Rao's Q (mean pairwise spectral distance)
  - Spectral CV (mean coefficient of variation across bands)
  - Spectral Shannon (mean entropy across binned bands)
  - PCA-based: FRic (convex hull volume), FDiv, FEve
  - Per-band mean and SD

Usage:
  python compute/compute_plot_spectral_1m.py
  python compute/compute_plot_spectral_1m.py --site HARV
  python compute/compute_plot_spectral_1m.py --grain 1   # only 1m grain
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import io
import re
import warnings
import zipfile

import numpy as np
import pandas as pd
import rasterio
from scipy.spatial import ConvexHull
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

from site_config import (
    SITES, VI_DIR, SPEC_DIV_DIR, VEG_STRUCT_DIR, CELL_SIZE,
)

warnings.filterwarnings("ignore")

PLOT_SIZE = 40  # meters
GRAIN_SIZES = [1, 2, 5, 10]  # meters
VI_BANDS = ["NDVI", "EVI", "ARVI", "PRI", "SAVI"]
ALL_BANDS = VI_BANDS + ["LAI", "fPAR"]
EXTRA_PRODUCTS = {
    "LAI":  ("DP3.30012.002", "_LAI.tif"),
    "fPAR": ("DP3.30014.002", "_fPAR.tif"),
}
MAX_SAMPLE_RAO = 500  # subsample for Rao's Q
MAX_SAMPLE_PCA = 2000  # subsample for PCA/ConvexHull/MST


# ─── Plot coordinate loading ────────────────────────────────────────────────

def load_plot_coords():
    """Load plot UTM coordinates from NEON vegetation structure data."""
    ppy_path = VEG_STRUCT_DIR / "vst_perplotperyear.csv"
    if not ppy_path.exists():
        print(f"ERROR: {ppy_path} not found")
        return pd.DataFrame()

    df = pd.read_csv(ppy_path, low_memory=False)

    # Need easting/northing columns
    coord_cols = [c for c in df.columns if c in
                  ["siteID", "plotID", "easting", "northing",
                   "decimalLatitude", "decimalLongitude"]]

    if "easting" in df.columns and "northing" in df.columns:
        coords = df[["siteID", "plotID", "easting", "northing"]].dropna()
        coords = coords.drop_duplicates(subset=["siteID", "plotID"])
        return coords

    print("WARNING: No UTM coordinates in vst_perplotperyear")
    return pd.DataFrame()


# ─── Tile discovery ──────────────────────────────────────────────────────────

def _parse_tile_coords(filename):
    m = re.search(r'_(\d{6})_(\d{7})_', filename)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None


def find_tile_for_point(site, x, y, year_tiles):
    """Find the tile key containing a given UTM point."""
    tile_e = int(x // 1000) * 1000
    tile_n = int(y // 1000) * 1000
    key = (tile_e, tile_n)
    return year_tiles.get(key)


def discover_tiles_by_year(site):
    """Discover VI tiles grouped by year, indexed by (easting, northing).

    Returns: {year: {(easting, northing): {"VI_zip": path, "LAI": path, ...}}}
    """
    vi_base = VI_DIR / "DP3.30026.002"
    result = {}

    for zp in vi_base.rglob(f"*_{site}_*_VegIndices.zip"):
        coords = _parse_tile_coords(zp.name)
        if not coords:
            continue
        # Extract year from path: .../FullSite/D01/2022_HARV_7/...
        m = re.search(r'[/\\](\d{4})[/\\]FullSite[/\\]', str(zp))
        if not m:
            continue
        year = int(m.group(1))

        if year not in result:
            result[year] = {}
        if coords not in result[year]:
            result[year][coords] = {}
        result[year][coords]["VI_zip"] = zp

    # Add LAI and fPAR tiles
    for bk, (dpid, suffix) in EXTRA_PRODUCTS.items():
        base = VI_DIR / dpid
        for tf in base.rglob(f"*_{site}_*{suffix}"):
            if "_error" in tf.name:
                continue
            coords = _parse_tile_coords(tf.name)
            if not coords:
                continue
            m = re.search(r'[/\\](\d{4})[/\\]FullSite[/\\]', str(tf))
            if not m:
                continue
            year = int(m.group(1))
            if year not in result:
                result[year] = {}
            if coords not in result[year]:
                result[year][coords] = {}
            result[year][coords][bk] = tf

    return result


# ─── 1m pixel extraction ────────────────────────────────────────────────────

def read_vi_band_from_zip(zip_path, band_name):
    """Read a single band from VI ZIP as (array, transform)."""
    with zipfile.ZipFile(str(zip_path)) as zf:
        target = [n for n in zf.namelist()
                  if n.endswith(f"_{band_name}.tif") and "_error" not in n]
        if not target:
            return None, None
        with zf.open(target[0]) as f:
            data = f.read()
        with rasterio.open(io.BytesIO(data)) as src:
            arr = src.read(1).astype(np.float32)
            if src.nodata is not None:
                arr[arr == src.nodata] = np.nan
            return arr, src.transform


def read_tif_band(path):
    with rasterio.open(str(path)) as src:
        arr = src.read(1).astype(np.float32)
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
        return arr, src.transform


def extract_plot_1m(tile_entry, plot_x, plot_y):
    """Extract all bands at 1m within a 40x40m plot footprint.

    Returns: (n_pixels, n_bands) array or None.
    """
    if "VI_zip" not in tile_entry:
        return None

    half = PLOT_SIZE / 2
    xmin = plot_x - half
    xmax = plot_x + half
    ymin = plot_y - half
    ymax = plot_y + half

    bands_data = []

    # VI bands from ZIP
    for band in VI_BANDS:
        arr, transform = read_vi_band_from_zip(tile_entry["VI_zip"], band)
        if arr is None:
            return None
        # Extract window
        row_start = int((transform.f - ymax) / abs(transform.e))
        row_end = int((transform.f - ymin) / abs(transform.e))
        col_start = int((xmin - transform.c) / transform.a)
        col_end = int((xmax - transform.c) / transform.a)

        # Clip to array bounds
        row_start = max(0, row_start)
        row_end = min(arr.shape[0], row_end)
        col_start = max(0, col_start)
        col_end = min(arr.shape[1], col_end)

        if row_end <= row_start or col_end <= col_start:
            return None

        patch = arr[row_start:row_end, col_start:col_end]
        bands_data.append(patch)

    # Extra bands (LAI, fPAR)
    for bk in ["LAI", "fPAR"]:
        if bk in tile_entry:
            arr, transform = read_tif_band(tile_entry[bk])
            row_start = max(0, int((transform.f - ymax) / abs(transform.e)))
            row_end = min(arr.shape[0], int((transform.f - ymin) / abs(transform.e)))
            col_start = max(0, int((xmin - transform.c) / transform.a))
            col_end = min(arr.shape[1], int((xmax - transform.c) / transform.a))
            if row_end > row_start and col_end > col_start:
                patch = arr[row_start:row_end, col_start:col_end]
            else:
                patch = np.full(bands_data[0].shape, np.nan, dtype=np.float32)
        else:
            patch = np.full(bands_data[0].shape, np.nan, dtype=np.float32)

        # Resize to match if needed
        target_shape = bands_data[0].shape
        if patch.shape != target_shape:
            patch = patch[:target_shape[0], :target_shape[1]]
            if patch.shape != target_shape:
                tmp = np.full(target_shape, np.nan, dtype=np.float32)
                h, w = min(patch.shape[0], target_shape[0]), min(patch.shape[1], target_shape[1])
                tmp[:h, :w] = patch[:h, :w]
                patch = tmp

        bands_data.append(patch)

    # Stack: (n_bands, H, W)
    stack = np.stack(bands_data, axis=0)
    return stack


# ─── Grain-size aggregation ─────────────────────────────────────────────────

def aggregate_to_grain(stack, grain_m):
    """Aggregate a (n_bands, H, W) 1m stack to a coarser grain.

    Returns (n_pixels, n_bands) 2D array.
    """
    if grain_m == 1:
        n_bands = stack.shape[0]
        flat = stack.reshape(n_bands, -1).T
        valid = ~np.any(np.isnan(flat), axis=1)
        return flat[valid]

    n_bands, h, w = stack.shape
    bh = (h // grain_m) * grain_m
    bw = (w // grain_m) * grain_m
    trimmed = stack[:, :bh, :bw]

    # Reshape to blocks
    reshaped = trimmed.reshape(n_bands, bh // grain_m, grain_m, bw // grain_m, grain_m)
    with np.errstate(all="ignore"):
        agg = np.nanmean(reshaped, axis=(2, 4))

    flat = agg.reshape(n_bands, -1).T
    valid = ~np.any(np.isnan(flat), axis=1)
    return flat[valid]


# ─── Diversity metrics ───────────────────────────────────────────────────────

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
    """PCA-based functional spectral diversity."""
    result = {"FRic": np.nan, "FDiv": np.nan, "FEve": np.nan}
    n = len(data)
    if n < n_pcs + 2:
        return result

    sample = data
    if n > MAX_SAMPLE_PCA:
        idx = np.random.choice(n, MAX_SAMPLE_PCA, replace=False)
        sample = data[idx]

    pca = PCA(n_components=min(n_pcs, sample.shape[1]))
    scores = pca.fit_transform(sample)

    # FRic: convex hull
    try:
        if scores.shape[1] >= 3 and len(scores) >= 4:
            hull = ConvexHull(scores[:, :3])
            result["FRic"] = float(hull.volume)
        elif scores.shape[1] >= 2 and len(scores) >= 3:
            hull = ConvexHull(scores[:, :2])
            result["FRic"] = float(hull.volume)
    except Exception:
        pass

    # FDiv
    centroid = np.mean(scores, axis=0)
    dists = np.sqrt(np.sum((scores - centroid) ** 2, axis=1))
    if np.max(dists) > 0:
        result["FDiv"] = float(np.mean(dists) / np.max(dists))

    # FEve: MST
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
            result["FEve"] = float(np.clip(eve, 0, 1))

    return result


def compute_all_metrics(data, grain_m):
    """Compute all spectral diversity metrics for a pixel array."""
    row = {
        "grain_m": grain_m,
        "n_pixels": len(data),
    }

    # Per-band stats
    for i, band in enumerate(ALL_BANDS):
        if i < data.shape[1]:
            row[f"{band}_mean"] = float(np.nanmean(data[:, i]))
            row[f"{band}_sd"] = float(np.nanstd(data[:, i]))

    # Diversity metrics
    row["rao_q"] = compute_rao_q(data)
    row["spectral_cv"] = compute_spectral_cv(data)
    row["spectral_shannon"] = compute_spectral_shannon(data)

    # PCA-based
    pca_metrics = compute_pca_diversity(data)
    row["spectral_FRic"] = pca_metrics["FRic"]
    row["spectral_FDiv"] = pca_metrics["FDiv"]
    row["spectral_FEve"] = pca_metrics["FEve"]

    return row


# ─── Main processing ────────────────────────────────────────────────────────

def process_site(site, plot_coords, grain_sizes):
    """Process all plots for a site across available years."""
    site_plots = plot_coords[plot_coords["siteID"] == site]
    if site_plots.empty:
        print(f"  SKIP {site}: no plot coordinates")
        return []

    tiles_by_year = discover_tiles_by_year(site)
    if not tiles_by_year:
        print(f"  SKIP {site}: no VI tiles")
        return []

    years = sorted(tiles_by_year.keys())
    print(f"  {site}: {len(site_plots)} plots, {len(years)} years")

    results = []
    for year in years:
        year_tiles = tiles_by_year[year]
        n_ok = 0

        for _, prow in site_plots.iterrows():
            px, py = prow["easting"], prow["northing"]
            tile_key = (int(px // 1000) * 1000, int(py // 1000) * 1000)
            tile_entry = year_tiles.get(tile_key)
            if tile_entry is None:
                continue

            # Extract 1m pixels
            stack = extract_plot_1m(tile_entry, px, py)
            if stack is None:
                continue

            # Compute at each grain size
            for grain in grain_sizes:
                data = aggregate_to_grain(stack, grain)
                if len(data) < 4:
                    continue

                row = compute_all_metrics(data, grain)
                row["siteID"] = site
                row["plotID"] = prow["plotID"]
                row["year"] = year
                results.append(row)

            n_ok += 1

        if n_ok > 0:
            print(f"    {year}: {n_ok}/{len(site_plots)} plots")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Plot-level spectral diversity at 1m resolution")
    parser.add_argument("--site", type=str, default=None)
    parser.add_argument("--grain", type=int, default=None,
                        help="Compute only this grain size (1, 2, 5, or 10)")
    args = parser.parse_args()

    SPEC_DIV_DIR.mkdir(parents=True, exist_ok=True)

    grain_sizes = [args.grain] if args.grain else GRAIN_SIZES

    # Load plot coordinates
    print("Loading plot coordinates...")
    plot_coords = load_plot_coords()
    if plot_coords.empty:
        print("ERROR: No plot coordinates found.")
        return
    print(f"  {len(plot_coords)} plots\n")

    sites = [args.site] if args.site else SITES
    all_results = []

    for site in sites:
        print(f"\n{'='*50}")
        rows = process_site(site, plot_coords, grain_sizes)
        all_results.extend(rows)

    if all_results:
        df = pd.DataFrame(all_results)
        out_path = SPEC_DIV_DIR / "plot_spectral_1m.csv"
        df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path} ({len(df)} rows)")

        # Summary by grain size
        print("\nPixels per plot by grain size:")
        for g in grain_sizes:
            sub = df[df["grain_m"] == g]
            if not sub.empty:
                print(f"  {g}m: {sub['n_pixels'].mean():.0f} pixels/plot, "
                      f"{len(sub)} plot-years, "
                      f"Rao's Q = {sub['rao_q'].mean():.4f}")


if __name__ == "__main__":
    main()
