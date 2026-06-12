"""
Spectral Diversity from NEON Vegetation Indices
=================================================
Computes spectral diversity metrics from downloaded NEON AOP vegetation
index products at 1m resolution, aggregated to 10m (Sentinel-2 aligned).

Input products (all 1m, 1km tiles):
  DP3.30026.002 (VI ZIPs):  NDVI, EVI, ARVI, PRI, SAVI
  DP3.30012.002 (LAI TIFs): LAI
  DP3.30014.002 (fPAR TIFs): fPAR

Output per site-year:
  - Multi-band 10m GeoTIFF: mean + SD of each index per 10m cell
  - Plot-level CSV: spectral diversity metrics (Rao's Q, CV, Shannon)
  - PCA-based functional spectral diversity (richness, divergence, evenness)

Usage:
  python compute_spectral_diversity.py                   # all sites
  python compute_spectral_diversity.py --site HARV       # single site
  python compute_spectral_diversity.py --dry-run         # list only
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import io
import re
import warnings
import zipfile
from collections import defaultdict

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA

from site_config import (
    SITES, SITE_EPSG, VI_DIR, SPEC_DIV_DIR, FSD_DIR, CELL_SIZE,
    get_fsd_files,
)

warnings.filterwarnings("ignore")

# ─── Band definitions ───────────────────────────────────────────────────────

VI_BANDS = ["NDVI", "EVI", "ARVI", "PRI", "SAVI"]  # from DP3.30026.002 ZIPs
EXTRA_BANDS = {
    "LAI": ("DP3.30012.002", "LAI", "_LAI.tif"),
    "fPAR": ("DP3.30014.002", "FPAR", "_fPAR.tif"),
}
ALL_BANDS = VI_BANDS + list(EXTRA_BANDS.keys())
N_BANDS_OUT = len(ALL_BANDS) * 2  # mean + SD per band


# ─── Tile discovery ─────────────────────────────────────────────────────────

def _parse_tile_coords(filename):
    """Extract (easting, northing) from NEON tile filename."""
    m = re.search(r'_(\d{6})_(\d{7})_', filename)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None


def discover_vi_zips(site):
    """Find all VI ZIP files for a site, grouped by (year, easting, northing)."""
    vi_base = VI_DIR / "DP3.30026.002"
    tiles = {}
    for zp in vi_base.rglob(f"*_{site}_*_VegIndices.zip"):
        coords = _parse_tile_coords(zp.name)
        if not coords:
            continue
        # Extract year from path
        m = re.search(r'(\d{4})', str(zp.parts[-5]) if len(zp.parts) > 5 else zp.name)
        if not m:
            continue
        year = int(m.group(1))
        tiles[(year, *coords)] = zp
    return tiles


def discover_extra_tifs(site, band_key):
    """Find TIF files for LAI or fPAR for a site."""
    dpid, subdir, suffix = EXTRA_BANDS[band_key]
    base = VI_DIR / dpid
    tiles = {}
    for tf in base.rglob(f"*_{site}_*{suffix}"):
        if "_error" in tf.name:
            continue
        coords = _parse_tile_coords(tf.name)
        if not coords:
            continue
        m = re.search(r'(\d{4})', str(tf.parts[-5]) if len(tf.parts) > 5 else tf.name)
        if not m:
            continue
        year = int(m.group(1))
        tiles[(year, *coords)] = tf
    return tiles


def discover_site_tiles(site):
    """Discover all available tiles for a site.

    Returns dict: {(year, easting, northing): {"VI_zip": path, "LAI": path, "fPAR": path}}
    """
    vi_zips = discover_vi_zips(site)

    extra = {}
    for bk in EXTRA_BANDS:
        extra[bk] = discover_extra_tifs(site, bk)

    # Merge by tile key
    all_keys = set(vi_zips.keys())
    for bk in extra:
        all_keys |= set(extra[bk].keys())

    tiles = {}
    for key in sorted(all_keys):
        entry = {}
        if key in vi_zips:
            entry["VI_zip"] = vi_zips[key]
        for bk in extra:
            if key in extra[bk]:
                entry[bk] = extra[bk][key]
        tiles[key] = entry

    return tiles


# ─── Tile reading ────────────────────────────────────────────────────────────

def read_vi_from_zip(zip_path, band_name):
    """Read a single VI band from a ZIP file as numpy array."""
    with zipfile.ZipFile(str(zip_path)) as zf:
        target = [n for n in zf.namelist()
                  if n.endswith(f"_{band_name}.tif") and "_error" not in n]
        if not target:
            return None, None
        with zf.open(target[0]) as f:
            data = f.read()
        with rasterio.open(io.BytesIO(data)) as src:
            arr = src.read(1).astype(np.float32)
            arr[arr == src.nodata] = np.nan
            return arr, src.transform


def read_tif(path):
    """Read a single-band TIF as numpy array."""
    with rasterio.open(str(path)) as src:
        arr = src.read(1).astype(np.float32)
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
        return arr, src.transform


def load_tile_stack(tile_entry):
    """Load all bands for a single tile into a (n_bands, H, W) stack.

    Returns (stack, transform) or (None, None) if missing VI data.
    """
    if "VI_zip" not in tile_entry:
        return None, None

    bands = []
    transform = None

    # VI bands from ZIP
    for band in VI_BANDS:
        arr, t = read_vi_from_zip(tile_entry["VI_zip"], band)
        if arr is None:
            return None, None
        bands.append(arr)
        if transform is None:
            transform = t

    # Extra bands (LAI, fPAR)
    for bk in EXTRA_BANDS:
        if bk in tile_entry:
            arr, _ = read_tif(tile_entry[bk])
            # Resize if needed (should be 1000x1000 like VI)
            if arr.shape != bands[0].shape:
                arr = arr[:bands[0].shape[0], :bands[0].shape[1]]
            bands.append(arr)
        else:
            bands.append(np.full(bands[0].shape, np.nan, dtype=np.float32))

    return np.stack(bands, axis=0), transform


# ─── 10m aggregation ────────────────────────────────────────────────────────

def aggregate_tile_10m(stack):
    """Aggregate a (n_bands, 1000, 1000) tile from 1m to 10m.

    Returns (mean_stack, sd_stack) each of shape (n_bands, 100, 100).
    """
    n_bands, h, w = stack.shape
    # Trim to nearest multiple of 10
    h10 = (h // 10) * 10
    w10 = (w // 10) * 10
    trimmed = stack[:, :h10, :w10]

    # Reshape to blocks of 10x10
    reshaped = trimmed.reshape(n_bands, h10 // 10, 10, w10 // 10, 10)

    with np.errstate(all="ignore"):
        mean_vals = np.nanmean(reshaped, axis=(2, 4))
        sd_vals = np.nanstd(reshaped, axis=(2, 4))

    return mean_vals, sd_vals


# ─── Spectral diversity metrics ─────────────────────────────────────────────

def rao_q_fast(stack_2d, n_sample=500):
    """Compute Rao's Q on a (n_pixels, n_bands) matrix.

    Uses random sampling for speed when n_pixels > n_sample.
    Q = sum(d_ij * p_i * p_j) where d = Euclidean distance.
    For uniform weights: Q = mean(d_ij) for all pairs.
    """
    valid = stack_2d[~np.any(np.isnan(stack_2d), axis=1)]
    n = len(valid)
    if n < 2:
        return np.nan

    if n > n_sample:
        idx = np.random.choice(n, n_sample, replace=False)
        valid = valid[idx]
        n = n_sample

    # Pairwise Euclidean distances (vectorized)
    diff = valid[:, np.newaxis, :] - valid[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diff ** 2, axis=2))
    # Mean of upper triangle
    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    return float(np.mean(dists[mask]))


def spectral_cv(stack_2d):
    """Mean CV across all bands."""
    valid = stack_2d[~np.any(np.isnan(stack_2d), axis=1)]
    if len(valid) < 2:
        return np.nan
    means = np.nanmean(valid, axis=0)
    sds = np.nanstd(valid, axis=0)
    cvs = np.where(np.abs(means) > 1e-6, sds / np.abs(means), np.nan)
    return float(np.nanmean(cvs))


def spectral_shannon(stack_2d, n_bins=20):
    """Shannon entropy across binned spectral values (mean across bands)."""
    valid = stack_2d[~np.any(np.isnan(stack_2d), axis=1)]
    if len(valid) < 2:
        return np.nan
    entropies = []
    for b in range(valid.shape[1]):
        vals = valid[:, b]
        counts, _ = np.histogram(vals, bins=n_bins)
        p = counts / counts.sum()
        p = p[p > 0]
        entropies.append(-np.sum(p * np.log(p)))
    return float(np.mean(entropies))


def functional_spectral_diversity(stack_2d, n_pcs=3):
    """PCA-based functional spectral diversity.

    Returns dict with FRic (convex hull volume), FDiv, FEve or NaN if insufficient data.
    """
    valid = stack_2d[~np.any(np.isnan(stack_2d), axis=1)]
    result = {"spectral_FRic": np.nan, "spectral_FDiv": np.nan, "spectral_FEve": np.nan}

    if len(valid) < n_pcs + 2:
        return result

    # Subsample for PCA-based metrics (avoid memory issues)
    MAX_PCA = 5000
    if len(valid) > MAX_PCA:
        idx = np.random.choice(len(valid), MAX_PCA, replace=False)
        valid = valid[idx]

    # PCA
    pca = PCA(n_components=min(n_pcs, valid.shape[1]))
    scores = pca.fit_transform(valid)

    # Functional richness: convex hull volume
    if scores.shape[1] >= 3 and len(scores) >= 4:
        try:
            hull = ConvexHull(scores[:, :3])
            result["spectral_FRic"] = float(hull.volume)
        except Exception:
            pass
    elif scores.shape[1] >= 2 and len(scores) >= 3:
        try:
            hull = ConvexHull(scores[:, :2])
            result["spectral_FRic"] = float(hull.volume)
        except Exception:
            pass

    # Functional divergence: mean distance to centroid / max distance
    centroid = np.mean(scores, axis=0)
    dists = np.sqrt(np.sum((scores - centroid) ** 2, axis=1))
    if np.max(dists) > 0:
        result["spectral_FDiv"] = float(np.mean(dists) / np.max(dists))

    # Functional evenness: regularity of MST edge lengths
    # Subsample to avoid memory explosion with pdist
    MAX_FEVE = 2000
    if len(scores) >= 3:
        from scipy.sparse.csgraph import minimum_spanning_tree
        from scipy.spatial.distance import pdist, squareform
        sc = scores
        if len(sc) > MAX_FEVE:
            idx = np.random.choice(len(sc), MAX_FEVE, replace=False)
            sc = sc[idx]
        dist_matrix = squareform(pdist(sc))
        mst = minimum_spanning_tree(dist_matrix)
        edges = mst.data[mst.data > 0]
        if len(edges) > 1:
            total = edges.sum()
            pev = edges / total
            n_edges = len(edges)
            expected = 1.0 / n_edges
            eve = 1.0 - np.sum(np.abs(pev - expected)) / (2.0 * (1.0 - expected))
            result["spectral_FEve"] = float(np.clip(eve, 0, 1))

    return result


# ─── Site processing ─────────────────────────────────────────────────────────

def process_site_year(site, year, tiles_for_year, dry_run=False):
    """Process all tiles for a site-year.

    Returns (10m_mosaic_bands, summary_dict) or None.
    """
    keys = sorted(tiles_for_year.keys())
    if not keys:
        return None

    if dry_run:
        print(f"    {year}: {len(keys)} tiles")
        return None

    # Determine mosaic extent from tile coordinates
    eastings = [k[1] for k in keys]
    northings = [k[2] for k in keys]
    xmin = min(eastings)
    ymin = min(northings)
    xmax = max(eastings) + 1000  # tiles are 1km
    ymax = max(northings) + 1000

    # 10m grid
    ncols_10m = (xmax - xmin) // CELL_SIZE
    nrows_10m = (ymax - ymin) // CELL_SIZE

    n_out = len(ALL_BANDS) * 2  # mean + SD
    mosaic = np.full((n_out, nrows_10m, ncols_10m), np.nan, dtype=np.float32)

    # Also collect all valid pixels for site-level PCA
    all_pixels = []
    n_tiles_ok = 0

    for key in keys:
        _, e, n = key
        entry = tiles_for_year[key]

        stack, transform = load_tile_stack(entry)
        if stack is None:
            continue

        mean_vals, sd_vals = aggregate_tile_10m(stack)
        n_tiles_ok += 1

        # Place into mosaic
        col_off = (e - xmin) // CELL_SIZE
        row_off = (ymax - n - 1000) // CELL_SIZE
        h, w = mean_vals.shape[1], mean_vals.shape[2]
        r_end = min(row_off + h, nrows_10m)
        c_end = min(col_off + w, ncols_10m)
        hh = r_end - row_off
        ww = c_end - col_off

        for b in range(len(ALL_BANDS)):
            mosaic[b * 2, row_off:r_end, col_off:c_end] = mean_vals[b, :hh, :ww]
            mosaic[b * 2 + 1, row_off:r_end, col_off:c_end] = sd_vals[b, :hh, :ww]

        # Sample pixels for PCA (subsample to avoid memory issues)
        flat = stack.reshape(stack.shape[0], -1).T  # (n_pixels, n_bands)
        valid_mask = ~np.any(np.isnan(flat), axis=1)
        valid_pixels = flat[valid_mask]
        if len(valid_pixels) > 1000:
            idx = np.random.choice(len(valid_pixels), 1000, replace=False)
            valid_pixels = valid_pixels[idx]
        all_pixels.append(valid_pixels)

    if n_tiles_ok == 0:
        return None

    # Site-level spectral diversity from sampled pixels
    if all_pixels:
        combined = np.vstack(all_pixels)
        site_rao = rao_q_fast(combined, n_sample=1000)
        site_cv = spectral_cv(combined)
        site_shannon = spectral_shannon(combined)
        site_func = functional_spectral_diversity(combined)
    else:
        site_rao = site_cv = site_shannon = np.nan
        site_func = {"spectral_FRic": np.nan, "spectral_FDiv": np.nan, "spectral_FEve": np.nan}

    summary = {
        "siteID": site,
        "year": year,
        "n_tiles": n_tiles_ok,
        "rao_q": site_rao,
        "spectral_cv": site_cv,
        "spectral_shannon": site_shannon,
        **site_func,
    }

    print(f"    {year}: {n_tiles_ok} tiles, "
          f"Rao's Q={site_rao:.4f}, Shannon={site_shannon:.3f}")

    return mosaic, xmin, ymax, ncols_10m, nrows_10m, summary


def process_site(site, dry_run=False):
    """Process all years for a site."""
    all_tiles = discover_site_tiles(site)
    if not all_tiles:
        print(f"  SKIP {site}: no VI tiles found")
        return []

    # Group by year
    by_year = defaultdict(dict)
    for (yr, e, n), entry in all_tiles.items():
        by_year[yr][(yr, e, n)] = entry

    years = sorted(by_year.keys())
    print(f"  {site}: {len(years)} years, {len(all_tiles)} total tiles")

    if dry_run:
        for yr in years:
            print(f"    {yr}: {len(by_year[yr])} tiles")
        return []

    SPEC_DIV_DIR.mkdir(parents=True, exist_ok=True)
    summaries = []

    for yr in years:
        result = process_site_year(site, yr, by_year[yr])
        if result is None:
            continue

        mosaic, xmin, ymax, ncols, nrows, summary = result
        summaries.append(summary)

        # Write 10m GeoTIFF
        epsg = SITE_EPSG.get(site)
        transform = from_origin(xmin, ymax, CELL_SIZE, CELL_SIZE)
        profile = {
            "driver": "GTiff",
            "dtype": "float32",
            "width": ncols,
            "height": nrows,
            "count": mosaic.shape[0],
            "transform": transform,
            "compress": "lzw",
            "tiled": True,
            "blockxsize": 256,
            "blockysize": 256,
            "nodata": np.nan,
        }
        if epsg:
            profile["crs"] = CRS.from_epsg(epsg)

        out_path = SPEC_DIV_DIR / f"{yr}_{site}_spectral_10m.tif"
        with rasterio.open(str(out_path), "w", **profile) as dst:
            for i in range(mosaic.shape[0]):
                dst.write(mosaic[i], i + 1)
                band_idx = i // 2
                stat = "mean" if i % 2 == 0 else "sd"
                dst.set_band_description(i + 1, f"{ALL_BANDS[band_idx]}_{stat}")

    return summaries


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compute spectral diversity from NEON VI data")
    parser.add_argument("--site", type=str, default=None, help="Process single site")
    parser.add_argument("--dry-run", action="store_true", help="List tiles only")
    args = parser.parse_args()

    sites = [args.site] if args.site else SITES

    print(f"Spectral Diversity Computation")
    print(f"Input: {VI_DIR}")
    print(f"Output: {SPEC_DIV_DIR}")
    print(f"Bands: {ALL_BANDS}")
    print(f"Sites: {len(sites)}\n")

    all_summaries = []
    for site in sites:
        print(f"\n{'='*50}")
        results = process_site(site, dry_run=args.dry_run)
        all_summaries.extend(results)

    # Save summary CSV
    if all_summaries and not args.dry_run:
        import pandas as pd
        df = pd.DataFrame(all_summaries)
        csv_path = SPEC_DIV_DIR / "spectral_diversity_all.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nSummary: {csv_path} ({len(df)} site-years)")


if __name__ == "__main__":
    main()
