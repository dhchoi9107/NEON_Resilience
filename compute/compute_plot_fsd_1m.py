"""
Plot-level FSD at 1m Resolution (40m plot + 50m buffer)
========================================================
Computes all 21 structural diversity metrics at 1m grid resolution
within each NEON plot footprint + 50m buffer (140x140m window).

Instead of computing 1m FSD across the entire site (extremely expensive),
this script only processes the area around each plot where the data is
actually needed for functional diversity analysis.

Output per plot: 140x140 = 19,600 pixels (vs 16 pixels at 10m)

Usage:
  python compute/compute_plot_fsd_1m.py
  python compute/compute_plot_fsd_1m.py --site HARV
  python compute/compute_plot_fsd_1m.py --dry-run
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import time
import warnings

import numpy as np
import pandas as pd
import laspy
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin

from site_config import (
    SITES, SITE_EPSG, NEON_BASE, VEG_STRUCT_DIR,
    FSD_BANDS,
)

warnings.filterwarnings("ignore")

# Import core metric functions from compute_fsd
# We set CELL_SIZE=1 for the imported module
import compute.compute_fsd as fsd_mod

PLOT_HALF = 20       # half of 40m plot
BUFFER = 50          # buffer around plot edge
WINDOW_HALF = PLOT_HALF + BUFFER  # 70m from plot center
WINDOW_SIZE = WINDOW_HALF * 2     # 140m total

CELL_SIZE_1M = 1
OUTPUT_DIR = Path("E:/neon_lidar/structural_diversity_1m_plots")


def load_plot_coords():
    """Load plot UTM coordinates."""
    ppy_path = VEG_STRUCT_DIR / "vst_perplotperyear.csv"
    if not ppy_path.exists():
        return pd.DataFrame()
    df = pd.read_csv(ppy_path, low_memory=False)
    if "easting" not in df.columns:
        return pd.DataFrame()
    coords = df[["siteID", "plotID", "easting", "northing"]].dropna()
    return coords.drop_duplicates(subset=["siteID", "plotID"])


def find_laz_files_for_plot(site, plot_x, plot_y):
    """Find LAZ tiles that overlap the plot + buffer window."""
    # Window bounds
    xmin = plot_x - WINDOW_HALF
    xmax = plot_x + WINDOW_HALF
    ymin = plot_y - WINDOW_HALF
    ymax = plot_y + WINDOW_HALF

    # Which 1km tiles could overlap?
    tile_xmin = int(xmin // 1000) * 1000
    tile_xmax = int(xmax // 1000) * 1000
    tile_ymin = int(ymin // 1000) * 1000
    tile_ymax = int(ymax // 1000) * 1000

    laz_files = []
    for pattern_base in [NEON_BASE / "neon-aop-products",
                         NEON_BASE / "neon-aop-provisional-products"]:
        for laz in pattern_base.rglob(f"*{site}*/ClassifiedPointCloud/*.laz"):
            # Parse tile coords from filename
            import re
            m = re.search(r'_(\d{6})_(\d{7})_', laz.name)
            if not m:
                continue
            te, tn = int(m.group(1)), int(m.group(2))
            # Check overlap
            if (te + 1000 > xmin and te < xmax and
                    tn + 1000 > ymin and tn < ymax):
                laz_files.append(laz)

    return laz_files


def process_plot_1m(laz_files, plot_x, plot_y):
    """Process LAZ files to compute 21-band FSD at 1m for a plot window.

    Returns (metrics_grid, xmin, ymax) or None.
    metrics_grid shape: (21, 140, 140)
    """
    # Window bounds
    xmin = plot_x - WINDOW_HALF
    xmax = plot_x + WINDOW_HALF
    ymin = plot_y - WINDOW_HALF
    ymax = plot_y + WINDOW_HALF

    # Read and merge points from all overlapping LAZ tiles
    all_x, all_y, all_z, all_cls = [], [], [], []

    for laz_path in laz_files:
        try:
            las = laspy.read(str(laz_path))
        except Exception:
            continue

        x = np.array(las.x, dtype=np.float64)
        y = np.array(las.y, dtype=np.float64)
        z = np.array(las.z, dtype=np.float64)
        cls = np.array(las.classification, dtype=np.uint8)

        # Spatial clip to window
        mask = (x >= xmin) & (x < xmax) & (y >= ymin) & (y < ymax)
        if mask.sum() == 0:
            continue

        all_x.append(x[mask])
        all_y.append(y[mask])
        all_z.append(z[mask])
        all_cls.append(cls[mask])

    if not all_x:
        return None

    x = np.concatenate(all_x)
    y = np.concatenate(all_y)
    z = np.concatenate(all_z)
    classification = np.concatenate(all_cls)

    if len(z) < 50:
        return None

    # Classification filter
    cls_mask = np.isin(classification, [1, 2, 3, 4, 5])
    x, y, z = x[cls_mask], y[cls_mask], z[cls_mask]
    classification = classification[cls_mask]

    if len(z) < 50:
        return None

    # Height normalization using compute_fsd's function
    z_norm = fsd_mod.normalize_heights(x, y, z, classification)
    if z_norm is None:
        return None

    valid = ~np.isnan(z_norm)
    x, y, z = x[valid], y[valid], z_norm[valid]

    # Height filter
    hm = (z >= -1.0) & (z <= 100.0)
    x, y, z = x[hm], y[hm], z[hm]

    if len(z) < 50:
        return None

    # Grid assignment at 1m
    n_cols = WINDOW_SIZE  # 140
    n_rows = WINDOW_SIZE  # 140
    n_cells = n_rows * n_cols

    col_idx = np.clip(((x - xmin) / CELL_SIZE_1M).astype(np.int32), 0, n_cols - 1)
    row_idx = np.clip((n_rows - 1 - ((y - ymin) / CELL_SIZE_1M).astype(np.int32)), 0, n_rows - 1)
    cell_id = (row_idx * n_cols + col_idx).astype(np.int32)

    # Outlier removal (per cell, 6-sigma)
    x_c, y_c, z_c, cid_c = fsd_mod._clean_and_assign(x, y, z, cell_id, n_cells)
    if len(z_c) == 0:
        return None

    # ── Compute all 21 metrics (replicating compute_fsd logic at 1m) ──
    Z0 = fsd_mod.Z0
    N_BANDS = fsd_mod.N_BANDS
    metrics = np.full((N_BANDS, n_rows, n_cols), np.nan, dtype=np.float32)

    count_all = np.bincount(cid_c, minlength=n_cells).astype(np.float64)
    sum_all = np.bincount(cid_c, weights=z_c, minlength=n_cells)
    sumsq_all = np.bincount(cid_c, weights=z_c**2, minlength=n_cells)

    has_pts = count_all > 0
    mean_all = np.zeros(n_cells)
    mean_all[has_pts] = sum_all[has_pts] / count_all[has_pts]
    var_all = np.zeros(n_cells)
    var_all[has_pts] = np.maximum(sumsq_all[has_pts] / count_all[has_pts] - mean_all[has_pts]**2, 0)
    sd_all = np.sqrt(var_all)

    # Vegetation subset
    veg_mask = z_c >= Z0
    z_veg = z_c[veg_mask]
    cid_veg = cid_c[veg_mask]

    count_veg = np.bincount(cid_veg, minlength=n_cells).astype(np.float64)
    sum_veg = np.bincount(cid_veg, weights=z_veg, minlength=n_cells)
    sumsq_veg = np.bincount(cid_veg, weights=z_veg**2, minlength=n_cells)

    has_veg = count_veg > 1
    mean_veg = np.full(n_cells, np.nan)
    mean_veg[has_veg] = sum_veg[has_veg] / count_veg[has_veg]
    var_veg = np.full(n_cells, np.nan)
    var_veg[has_veg] = np.maximum(
        (sumsq_veg[has_veg] - sum_veg[has_veg]**2 / count_veg[has_veg]) / (count_veg[has_veg] - 1), 0)
    sd_veg = np.sqrt(var_veg)

    # Band 0: rumple
    r_val = np.full(n_cells, np.nan)
    m_safe = has_pts & (mean_all > 0.01)
    r_val[m_safe] = 1.0 + sd_all[m_safe] / mean_all[m_safe]
    metrics[0] = r_val.reshape(n_rows, n_cols)

    # Band 1: top_rugosity
    z_rug = z_c.copy()
    z_rug[z_rug < Z0] = 0
    sum_rug = np.bincount(cid_c, weights=z_rug, minlength=n_cells)
    sumsq_rug = np.bincount(cid_c, weights=z_rug**2, minlength=n_cells)
    mean_rug = np.zeros(n_cells)
    mean_rug[has_pts] = sum_rug[has_pts] / count_all[has_pts]
    var_rug = np.zeros(n_cells)
    var_rug[has_pts] = np.maximum(sumsq_rug[has_pts] / count_all[has_pts] - mean_rug[has_pts]**2, 0)
    sd_rug = np.full(n_cells, np.nan)
    sd_rug[has_pts] = np.sqrt(var_rug[has_pts])
    metrics[1] = sd_rug.reshape(n_rows, n_cols)

    # Bands 2-7: mean_max_canopy_ht, max_canopy_ht, deepgap, meanH, vert_sd, vertCV
    metrics[2] = mean_veg.reshape(n_rows, n_cols)
    order_max = np.lexsort((z_c, cid_c))
    z_sm = z_c[order_max]
    cid_sm = cid_c[order_max]
    _, last_idx = np.unique(cid_sm[::-1], return_index=True)
    last_idx = len(z_sm) - 1 - last_idx
    max_ht = np.full(n_cells, np.nan, dtype=np.float32)
    max_ht[cid_sm[last_idx]] = z_sm[last_idx]
    metrics[3] = max_ht.reshape(n_rows, n_cols)

    count_ground = np.bincount(cid_c, weights=(z_rug == 0).astype(np.float64), minlength=n_cells)
    dgf = np.full(n_cells, np.nan)
    dgf[has_pts] = count_ground[has_pts] / count_all[has_pts]
    metrics[4] = dgf.reshape(n_rows, n_cols)
    metrics[5] = mean_veg.reshape(n_rows, n_cols)
    metrics[6] = sd_veg.reshape(n_rows, n_cols)
    vertCV = np.full(n_cells, np.nan)
    cv_mask = has_veg & (mean_veg > 0)
    vertCV[cv_mask] = sd_veg[cv_mask] / mean_veg[cv_mask]
    metrics[7] = vertCV.reshape(n_rows, n_cols)

    # Band 10: Gini
    metrics[10] = fsd_mod._per_cell_sorted_metric(
        z_veg, cid_veg, n_cells, fsd_mod._gini).reshape(n_rows, n_cols)

    # Bands 11-12: GFP, VCI
    gfp_arr = np.full(n_cells, np.nan, dtype=np.float32)
    vci_arr = np.full(n_cells, np.nan, dtype=np.float32)
    order = np.argsort(cid_c)
    z_s = z_c[order]
    cid_s = cid_c[order]
    if len(cid_s) > 0:
        changes = np.where(np.diff(cid_s) != 0)[0] + 1
        starts = np.concatenate([[0], changes])
        ends = np.concatenate([changes, [len(z_s)]])
        gids = cid_s[starts]
        for i in range(len(starts)):
            s, e = starts[i], ends[i]
            if e - s >= 2:
                g, v = fsd_mod._gfp_vci(z_s[s:e])
                gfp_arr[gids[i]] = g
                vci_arr[gids[i]] = v
    metrics[11] = gfp_arr.reshape(n_rows, n_cols)
    metrics[12] = vci_arr.reshape(n_rows, n_cols)

    # Bands 13-16: quantiles
    if len(z_veg) > 0:
        for qi, q in enumerate([25, 50, 75, 95]):
            def _qfunc(zg, _q=q):
                return np.percentile(zg, _q)
            metrics[13 + qi] = fsd_mod._per_cell_sorted_metric(
                z_veg, cid_veg, n_cells, _qfunc).reshape(n_rows, n_cols)

    # Band 17: HeightRatio
    count_b10 = np.bincount(cid_c, weights=(z_c < 10).astype(np.float64), minlength=n_cells)
    hr = np.full(n_cells, np.nan)
    hr[has_pts] = count_b10[has_pts] / count_all[has_pts] * 100
    metrics[17] = hr.reshape(n_rows, n_cols)

    # Bands 18-20: FHD, LAI, LAI_subcanopy
    fhd_arr = np.full(n_cells, np.nan, dtype=np.float32)
    lai_arr = np.full(n_cells, np.nan, dtype=np.float32)
    laisub_arr = np.full(n_cells, np.nan, dtype=np.float32)
    order2 = np.argsort(cid_c)
    z_s2 = z_c[order2]
    cid_s2 = cid_c[order2]
    if len(cid_s2) > 0:
        changes2 = np.where(np.diff(cid_s2) != 0)[0] + 1
        starts2 = np.concatenate([[0], changes2])
        ends2 = np.concatenate([changes2, [len(z_s2)]])
        gids2 = cid_s2[starts2]
        for i in range(len(starts2)):
            s, e = starts2[i], ends2[i]
            if e - s >= 5:
                fhd_v, lai_v, laisub_v = fsd_mod._fhd_lai_lai_sub(z_s2[s:e])
                fhd_arr[gids2[i]] = fhd_v
                lai_arr[gids2[i]] = lai_v
                laisub_arr[gids2[i]] = laisub_v
    metrics[18] = fhd_arr.reshape(n_rows, n_cols)
    metrics[19] = lai_arr.reshape(n_rows, n_cols)
    metrics[20] = laisub_arr.reshape(n_rows, n_cols)

    # Bands 8-9: mean_sd, sd_sd (skip sub-grid at 1m — not meaningful)
    # At 1m, sub-grid doesn't apply. Set to NaN.
    metrics[8] = np.full((n_rows, n_cols), np.nan, dtype=np.float32)
    metrics[9] = np.full((n_rows, n_cols), np.nan, dtype=np.float32)

    return metrics, xmin, ymax


def process_site(site, plot_coords):
    """Process all plots for a site, finding LAZ files per year."""
    site_plots = plot_coords[plot_coords["siteID"] == site]
    if site_plots.empty:
        print(f"  SKIP {site}: no plot coordinates")
        return

    # Discover available years for this site
    import re
    year_laz = {}  # {year: [laz_paths]}
    for laz in NEON_BASE.rglob(f"*{site}*classified*.laz"):
        if laz.name.startswith("._"):
            continue
        m = re.search(r'[/\\](\d{4})[/\\]FullSite[/\\]', str(laz))
        if m:
            yr = int(m.group(1))
            if yr not in year_laz:
                year_laz[yr] = []
            year_laz[yr].append(laz)

    if not year_laz:
        print(f"  SKIP {site}: no LAZ files")
        return

    years = sorted(year_laz.keys())
    print(f"  {site}: {len(site_plots)} plots, {len(years)} years")

    epsg = SITE_EPSG.get(site)
    site_dir = OUTPUT_DIR / site
    site_dir.mkdir(parents=True, exist_ok=True)

    for yr in years:
        t0 = time.time()
        n_ok = 0

        for _, prow in site_plots.iterrows():
            px, py = prow["easting"], prow["northing"]
            plot_id = prow["plotID"]

            out_path = site_dir / f"{yr}_{plot_id}_FSD_1m.tif"
            if out_path.exists():
                n_ok += 1
                continue

            # Find overlapping LAZ files for this plot
            xmin = px - WINDOW_HALF
            xmax = px + WINDOW_HALF
            ymin_p = py - WINDOW_HALF
            ymax_p = py + WINDOW_HALF

            plot_laz = []
            for laz in year_laz[yr]:
                m = re.search(r'_(\d{6})_(\d{7})_', laz.name)
                if not m:
                    continue
                te, tn = int(m.group(1)), int(m.group(2))
                if (te + 1000 > xmin and te < xmax and
                        tn + 1000 > ymin_p and tn < ymax_p):
                    plot_laz.append(laz)

            if not plot_laz:
                continue

            result = process_plot_1m(plot_laz, px, py)
            if result is None:
                continue

            grid, grid_xmin, grid_ymax = result

            # Save as small GeoTIFF
            transform = from_origin(grid_xmin, grid_ymax, CELL_SIZE_1M, CELL_SIZE_1M)
            profile = {
                "driver": "GTiff",
                "dtype": "float32",
                "width": WINDOW_SIZE,
                "height": WINDOW_SIZE,
                "count": len(FSD_BANDS),
                "transform": transform,
                "compress": "lzw",
                "nodata": np.nan,
            }
            if epsg:
                profile["crs"] = CRS.from_epsg(epsg)

            with rasterio.open(str(out_path), "w", **profile) as dst:
                for i in range(len(FSD_BANDS)):
                    dst.write(grid[i], i + 1)
                    dst.set_band_description(i + 1, FSD_BANDS[i])

            n_ok += 1

        elapsed = time.time() - t0
        print(f"    {yr}: {n_ok}/{len(site_plots)} plots ({elapsed:.0f}s)")


def main():
    parser = argparse.ArgumentParser(
        description="Compute plot-level FSD at 1m (40m plot + 50m buffer)")
    parser.add_argument("--site", type=str, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading plot coordinates...")
    plot_coords = load_plot_coords()
    if plot_coords.empty:
        print("ERROR: No plot coordinates. Run neon_veg_structure_download.py first.")
        return
    print(f"  {len(plot_coords)} plots\n")

    sites = [args.site] if args.site else SITES

    if args.dry_run:
        for site in sites:
            sp = plot_coords[plot_coords["siteID"] == site]
            print(f"  {site}: {len(sp)} plots")
        return

    for site in sites:
        print(f"\n{'='*50}")
        process_site(site, plot_coords)

    print(f"\nDone. Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
