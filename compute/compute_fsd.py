"""
NEON LiDAR Structural Diversity Metrics Calculator (v3)
=======================================================
- Height normalization via ground-point TIN interpolation (Class 2)
- 10m grid aligned to Sentinel-2 UTM grid (coordinates snapped to 10m multiples)
- 21-band multi-band GeoTIFF output per site-year

Metrics (21 bands):
  CHM:        rumple, top_rugosity, mean_max_canopy_ht, max_canopy_ht, deepgap_fraction
  Structural: meanH, vert_sd, vertCV, mean_sd, sd_sd, GC (Gini coefficient)
  LAI:        GFP, VCI, FHD, LAI, LAI_subcanopy
  Quantiles:  q25, q50, q75, q95
  Additional: HeightRatio

Usage:
  python compute_fsd.py                    # process all site-years
  python compute_fsd.py --site HARV        # process only HARV
  python compute_fsd.py --year 2022        # process only 2022
  python compute_fsd.py --workers 4        # set parallel workers
"""

import re
import argparse
import time
import warnings
from pathlib import Path

import numpy as np
import laspy
import rasterio
from rasterio.transform import from_origin
from rasterio.crs import CRS
from scipy.interpolate import LinearNDInterpolator
from joblib import Parallel, delayed
from tqdm import tqdm

warnings.filterwarnings("ignore")

# ─── Configuration ───────────────────────────────────────────────────────────
NEON_BASE = Path("E:/neon_lidar/DP1.30003.001")
OUTPUT_DIR = Path("E:/neon_lidar/structural_diversity")
CELL_SIZE = 10       # Sentinel-2 aligned
Z0 = 0.5             # minimum vegetation height cutoff
SIGMA_MULT = 6       # outlier removal: mean + 6*sd
SUBCELL_SIZE = 3     # 3m sub-grid for sd_sd (matches R code's 9m² grid)
GROUND_THIN = 5      # subsample ground points to every Nth for TIN speed

BEER_LAMBERT_K = 1.0  # extinction coefficient for LAD calculation
LAI_SUB_MAX = 5.0     # subcanopy max height for LAI_subcanopy

BAND_NAMES = [
    "rumple", "top_rugosity", "mean_max_canopy_ht", "max_canopy_ht", "deepgap_fraction",
    "meanH", "vert_sd", "vertCV", "mean_sd", "sd_sd", "GC",
    "GFP", "VCI",
    "q25", "q50", "q75", "q95",
    "HeightRatio",
    "FHD", "LAI", "LAI_subcanopy",
]
N_BANDS = len(BAND_NAMES)


# ─── Height Normalization ────────────────────────────────────────────────────

def normalize_heights(x, y, z, classification):
    """Normalize Z to height above ground using Class-2 ground points.
    Uses TIN (Delaunay triangulation + linear interpolation).
    Returns normalized Z array. Points below ground or without ground reference → NaN.
    """
    ground_mask = classification == 2
    n_ground = np.sum(ground_mask)

    if n_ground < 3:
        return None  # not enough ground points

    xg = x[ground_mask]
    yg = y[ground_mask]
    zg = z[ground_mask]

    # Thin ground points for speed (keep every Nth point, shuffled)
    if n_ground > 50000:
        step = max(1, n_ground // 50000)
        idx = np.arange(0, n_ground, step)
        xg, yg, zg = xg[idx], yg[idx], zg[idx]

    # Build TIN interpolator
    try:
        interp = LinearNDInterpolator(np.column_stack([xg, yg]), zg)
    except Exception:
        return None

    # Interpolate ground elevation at all points
    ground_elev = interp(np.column_stack([x, y]))

    # Normalize
    z_norm = z - ground_elev

    return z_norm


# ─── Vectorized Metric Helpers ───────────────────────────────────────────────

def _clean_and_assign(x, y, z, cell_ids, n_cells):
    """Filter z >= 0, apply per-cell 6σ outlier removal. Returns x, y, z, cell_ids."""
    mask0 = z >= 0
    z0 = z[mask0]
    cid0 = cell_ids[mask0]

    count = np.bincount(cid0, minlength=n_cells).astype(np.float64)
    sumz = np.bincount(cid0, weights=z0, minlength=n_cells)
    sumsq = np.bincount(cid0, weights=z0**2, minlength=n_cells)

    valid = count > 0
    mean = np.zeros(n_cells)
    mean[valid] = sumz[valid] / count[valid]
    var = np.zeros(n_cells)
    var[valid] = sumsq[valid] / count[valid] - mean[valid]**2
    var = np.maximum(var, 0.0)
    sd = np.sqrt(var)
    threshold = mean + SIGMA_MULT * sd

    keep = mask0.copy()
    keep[mask0] &= z0 < threshold[cid0]
    return x[keep], y[keep], z[keep], cell_ids[keep]


def _per_cell_sorted_metric(z, cid, n_cells, func):
    """Apply func(sorted_z_group) to each cell group. Returns (n_cells,) float32 array."""
    result = np.full(n_cells, np.nan, dtype=np.float32)
    order = np.lexsort((z, cid))
    z_s = z[order]
    cid_s = cid[order]

    if len(cid_s) == 0:
        return result

    changes = np.where(np.diff(cid_s) != 0)[0] + 1
    starts = np.concatenate([[0], changes])
    ends = np.concatenate([changes, [len(z_s)]])
    group_ids = cid_s[starts]

    for i in range(len(starts)):
        s, e = starts[i], ends[i]
        if e - s >= 2:
            result[group_ids[i]] = func(z_s[s:e])
    return result


def _gini(z_sorted):
    n = len(z_sorted)
    idx = np.arange(1, n + 1)
    sx = np.sum(z_sorted)
    if sx <= 0:
        return np.nan
    return (2 * np.sum(z_sorted * idx) / sx - (n + 1)) / (n - 1)


def _gfp_vci(z_group):
    """Returns (gfp, vci) tuple for one cell."""
    zv = z_group[z_group >= Z0]
    n = len(zv)
    if n < 2:
        return np.nan, np.nan
    zmax = np.max(zv)

    # GFP
    bins = np.arange(Z0, zmax + 1.0, 1.0)
    gfp_val = np.nan
    if len(bins) >= 2:
        counts, _ = np.histogram(zv, bins=bins)
        cum = 0
        gf_vals = []
        for c in counts:
            n_above = n - cum
            if n_above > 0:
                gf_vals.append(1.0 - c / n_above)
            cum += c
        if gf_vals:
            gfp_val = np.mean(gf_vals)

    # VCI
    vci_val = np.nan
    bins_v = np.arange(0, min(zmax, 100.0) + 1.0, 1.0)
    if len(bins_v) >= 2:
        counts_v, _ = np.histogram(zv, bins=bins_v)
        nv = np.sum(counts_v)
        if nv > 0:
            p = counts_v[counts_v > 0] / nv
            H = -np.sum(p * np.log(p))
            Hmax = np.log(len(counts_v))
            if Hmax > 0:
                vci_val = H / Hmax

    return gfp_val, vci_val


def _fhd_lai_lai_sub(z_group):
    """Compute FHD, LAI, LAI_subcanopy for one cell using Beer-Lambert LAD.
    Matches R code's lad.voxels2 + leafR::FHD + leafR::lai logic.
    Returns (fhd, lai, lai_sub) tuple.
    """
    z = z_group[z_group >= 0]
    if len(z) < 5:
        return np.nan, np.nan, np.nan

    zmax = np.max(z)
    if zmax < Z0:
        return np.nan, np.nan, np.nan

    maxZ = int(np.floor(zmax))
    if maxZ < 1:
        return np.nan, np.nan, np.nan

    # Count points per 1m height bin (bottom to top)
    bins = np.arange(0, maxZ + 1, 1)
    if len(bins) < 2:
        return np.nan, np.nan, np.nan
    counts, _ = np.histogram(z, bins=bins)
    n_layers = len(counts)

    if n_layers < 2:
        return np.nan, np.nan, np.nan

    # Reverse order: top to bottom (matches R code)
    counts_rev = counts[::-1].astype(np.float64)
    total = np.sum(counts_rev)
    if total == 0:
        return np.nan, np.nan, np.nan

    # Cumulative sum from top
    cumsum = np.cumsum(counts_rev)

    # pulse_out = total - cumsum
    pulse_out = total - cumsum

    # pulse_in: for first layer = total, for subsequent = pulse_out of previous
    pulse_in = np.concatenate([[total], pulse_out[:-1]])

    # LAD = log(pulse_in / pulse_out) * (1/k) * (1/dz)
    # Avoid division by zero or log of non-positive
    valid = (pulse_in > 0) & (pulse_out > 0)
    lad = np.full(n_layers, np.nan)
    lad[valid] = np.log(pulse_in[valid] / pulse_out[valid]) / BEER_LAMBERT_K

    # Remove last layer (unreliable, matches R code)
    lad = lad[:-1]
    if len(lad) == 0:
        return np.nan, np.nan, np.nan

    # Reverse back to bottom-to-top order, height labels
    lad = lad[::-1]
    heights = np.arange(0, len(lad))  # 0, 1, 2, ... m

    # Remove NaN and negative LAD
    lad_clean = lad.copy()
    lad_clean[np.isnan(lad_clean)] = 0
    lad_clean[lad_clean < 0] = 0

    # ── LAI: sum of LAD above z0 ──
    z0_idx = max(0, int(np.ceil(Z0)))
    lai = np.sum(lad_clean[z0_idx:])

    # ── LAI_subcanopy: sum of LAD from z0 to LAI_SUB_MAX ──
    sub_max_idx = min(len(lad_clean), int(LAI_SUB_MAX))
    lai_sub = np.sum(lad_clean[z0_idx:sub_max_idx])

    # ── FHD: Shannon entropy of relative LAD profile above z0 ──
    lad_above = lad_clean[z0_idx:]
    total_lad = np.sum(lad_above)
    if total_lad <= 0 or len(lad_above) < 2:
        fhd = np.nan
    else:
        p = lad_above / total_lad
        p = p[p > 0]
        fhd = -np.sum(p * np.log(p))

    return fhd, lai, lai_sub


# ─── Tile Processing ─────────────────────────────────────────────────────────

def parse_tile_coords(filename):
    m = re.search(r'_(\d{6})_(\d{7})_', filename)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def process_tile(laz_path):
    """Process one LAZ tile → (metrics_grid, epsg, bounds) or None."""
    try:
        las = laspy.read(str(laz_path))
    except Exception:
        return None

    x = np.array(las.x, dtype=np.float64)
    y = np.array(las.y, dtype=np.float64)
    z = np.array(las.z, dtype=np.float64)
    classification = np.array(las.classification, dtype=np.uint8)

    if len(z) < 10:
        return None

    # ── CRS ──
    epsg = None
    try:
        crs_obj = las.header.parse_crs()
        if crs_obj is not None:
            epsg = crs_obj.to_epsg()
    except Exception:
        pass
    # Fallback: infer UTM EPSG from coordinates (NEON sites are all in NAD83/WGS84 UTM)
    if epsg is None:
        try:
            lon_approx = np.median(x)
            lat_approx = np.median(y)
            # UTM coordinates: easting ~100k-900k, northing ~0-10M
            # If x looks like longitude (-180 to 180), compute UTM zone
            if -180 <= lon_approx <= 180:
                zone = int((lon_approx + 180) / 6) + 1
                epsg = 32600 + zone if lat_approx >= 0 else 32700 + zone
            elif 100_000 < lon_approx < 900_000:
                # Already in UTM — try WKT from VLR records
                for vlr in las.header.vlrs:
                    try:
                        txt = vlr.record_data.decode("ascii", errors="ignore")
                        import re as _re
                        m = _re.search(r'UTM\s+[Zz]one\s+(\d+)([NS])?', txt)
                        if m:
                            zone = int(m.group(1))
                            hemisphere = m.group(2)
                            epsg = 32600 + zone if hemisphere != 'S' else 32700 + zone
                            break
                    except Exception:
                        continue
        except Exception:
            pass
    if epsg is None:
        print(f"    WARNING: Could not determine CRS for {laz_path.name}")

    # ── Classification filter: keep only valid classes (1-5) ──
    # 1=unclassified, 2=ground, 3=low veg, 4=med veg, 5=high veg
    # Removes: 6=building, 7=noise, 18=high noise, etc.
    cls_mask = np.isin(classification, [1, 2, 3, 4, 5])
    x = x[cls_mask]
    y = y[cls_mask]
    z = z[cls_mask]
    classification = classification[cls_mask]

    if len(z) < 10:
        return None

    # ── Height normalization ──
    z_norm = normalize_heights(x, y, z, classification)
    if z_norm is None:
        return None

    # Remove points with NaN height (outside ground TIN convex hull)
    valid_mask = ~np.isnan(z_norm)
    x = x[valid_mask]
    y = y[valid_mask]
    z = z_norm[valid_mask]

    if len(z) < 10:
        return None

    # ── Abnormal height filter: remove Z < -1m or Z > 100m ──
    height_mask = (z >= -1.0) & (z <= 100.0)
    x = x[height_mask]
    y = y[height_mask]
    z = z[height_mask]

    if len(z) < 10:
        return None

    # ── Tile bounds snapped to Sentinel-2 grid (10m multiples) ──
    easting, northing = parse_tile_coords(laz_path.name)
    if easting is None:
        x_min = np.floor(np.min(x) / CELL_SIZE) * CELL_SIZE
        x_max = np.ceil(np.max(x) / CELL_SIZE) * CELL_SIZE
        y_min = np.floor(np.min(y) / CELL_SIZE) * CELL_SIZE
        y_max = np.ceil(np.max(y) / CELL_SIZE) * CELL_SIZE
    else:
        # NEON tiles are 1km, already 10m-aligned (multiples of 1000)
        x_min, x_max = float(easting), float(easting + 1000)
        y_min, y_max = float(northing), float(northing + 1000)

    n_cols = int((x_max - x_min) / CELL_SIZE)
    n_rows = int((y_max - y_min) / CELL_SIZE)
    if n_cols == 0 or n_rows == 0:
        return None

    n_cells = n_rows * n_cols

    # ── Assign points to grid cells ──
    col_idx = np.clip(((x - x_min) / CELL_SIZE).astype(np.int32), 0, n_cols - 1)
    row_idx = np.clip((n_rows - 1 - ((y - y_min) / CELL_SIZE).astype(np.int32)), 0, n_rows - 1)
    cell_id = (row_idx * n_cols + col_idx).astype(np.int32)

    # ── Clean (z>=0, outlier removal) ──
    x_clean, y_clean, z_clean, cid_clean = _clean_and_assign(x, y, z, cell_id, n_cells)
    if len(z_clean) == 0:
        return None

    metrics = np.full((N_BANDS, n_rows, n_cols), np.nan, dtype=np.float32)

    # ── Bulk bincount stats ──
    count_all = np.bincount(cid_clean, minlength=n_cells).astype(np.float64)
    sum_all = np.bincount(cid_clean, weights=z_clean, minlength=n_cells)
    sumsq_all = np.bincount(cid_clean, weights=z_clean**2, minlength=n_cells)

    has_pts = count_all > 0
    mean_all = np.zeros(n_cells)
    mean_all[has_pts] = sum_all[has_pts] / count_all[has_pts]
    var_all = np.zeros(n_cells)
    var_all[has_pts] = np.maximum(sumsq_all[has_pts] / count_all[has_pts] - mean_all[has_pts]**2, 0)
    sd_all = np.sqrt(var_all)

    # ── Vegetation subset (z >= z0) ──
    veg_mask = z_clean >= Z0
    z_veg = z_clean[veg_mask]
    cid_veg = cid_clean[veg_mask]

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

    # Band 0: rumple (distribution-based approximation)
    r_val = np.full(n_cells, np.nan)
    m_safe = has_pts & (mean_all > 0.01)
    r_val[m_safe] = 1.0 + sd_all[m_safe] / mean_all[m_safe]
    metrics[0] = r_val.reshape(n_rows, n_cols).astype(np.float32)

    # Band 1: top_rugosity
    z_rug = z_clean.copy()
    z_rug[z_rug < Z0] = 0
    sum_rug = np.bincount(cid_clean, weights=z_rug, minlength=n_cells)
    sumsq_rug = np.bincount(cid_clean, weights=z_rug**2, minlength=n_cells)
    mean_rug = np.zeros(n_cells)
    mean_rug[has_pts] = sum_rug[has_pts] / count_all[has_pts]
    var_rug = np.zeros(n_cells)
    var_rug[has_pts] = np.maximum(sumsq_rug[has_pts] / count_all[has_pts] - mean_rug[has_pts]**2, 0)
    sd_rug = np.full(n_cells, np.nan)
    sd_rug[has_pts] = np.sqrt(var_rug[has_pts])
    metrics[1] = sd_rug.reshape(n_rows, n_cols).astype(np.float32)

    # Band 2: mean_max_canopy_ht
    metrics[2] = mean_veg.reshape(n_rows, n_cols).astype(np.float32)

    # Band 3: max_canopy_ht (per-cell max via sorted unique)
    order_max = np.lexsort((z_clean, cid_clean))
    z_sm = z_clean[order_max]
    cid_sm = cid_clean[order_max]
    _, last_idx = np.unique(cid_sm[::-1], return_index=True)
    last_idx = len(z_sm) - 1 - last_idx
    max_ht = np.full(n_cells, np.nan, dtype=np.float32)
    max_ht[cid_sm[last_idx]] = z_sm[last_idx]
    metrics[3] = max_ht.reshape(n_rows, n_cols)

    # Band 4: deepgap_fraction
    count_ground = np.bincount(cid_clean, weights=(z_rug == 0).astype(np.float64), minlength=n_cells)
    dgf = np.full(n_cells, np.nan)
    dgf[has_pts] = count_ground[has_pts] / count_all[has_pts]
    metrics[4] = dgf.reshape(n_rows, n_cols).astype(np.float32)

    # Band 5: meanH
    metrics[5] = mean_veg.reshape(n_rows, n_cols).astype(np.float32)

    # Band 6: vert_sd
    metrics[6] = sd_veg.reshape(n_rows, n_cols).astype(np.float32)

    # Band 7: vertCV
    vertCV = np.full(n_cells, np.nan)
    cv_mask = has_veg & (mean_veg > 0)
    vertCV[cv_mask] = sd_veg[cv_mask] / mean_veg[cv_mask]
    metrics[7] = vertCV.reshape(n_rows, n_cols).astype(np.float32)

    # Band 10: Gini coefficient
    metrics[10] = _per_cell_sorted_metric(z_veg, cid_veg, n_cells, _gini).reshape(n_rows, n_cols)

    # Bands 11-12: GFP, VCI
    gfp_arr = np.full(n_cells, np.nan, dtype=np.float32)
    vci_arr = np.full(n_cells, np.nan, dtype=np.float32)
    order = np.argsort(cid_clean)
    z_s = z_clean[order]
    cid_s = cid_clean[order]
    if len(cid_s) > 0:
        changes = np.where(np.diff(cid_s) != 0)[0] + 1
        starts = np.concatenate([[0], changes])
        ends = np.concatenate([changes, [len(z_s)]])
        gids = cid_s[starts]
        for i in range(len(starts)):
            s, e = starts[i], ends[i]
            if e - s >= 2:
                g, v = _gfp_vci(z_s[s:e])
                gfp_arr[gids[i]] = g
                vci_arr[gids[i]] = v
    metrics[11] = gfp_arr.reshape(n_rows, n_cols)
    metrics[12] = vci_arr.reshape(n_rows, n_cols)

    # Bands 13-16: quantiles
    if len(z_veg) > 0:
        for qi, q in enumerate([25, 50, 75, 95]):
            def _qfunc(zg, _q=q):
                return np.percentile(zg, _q)
            metrics[13 + qi] = _per_cell_sorted_metric(z_veg, cid_veg, n_cells, _qfunc).reshape(n_rows, n_cols)

    # Band 17: HeightRatio
    count_b10 = np.bincount(cid_clean, weights=(z_clean < 10).astype(np.float64), minlength=n_cells)
    hr = np.full(n_cells, np.nan)
    hr[has_pts] = count_b10[has_pts] / count_all[has_pts] * 100
    metrics[17] = hr.reshape(n_rows, n_cols).astype(np.float32)

    # Bands 18-20: FHD, LAI, LAI_subcanopy (Beer-Lambert LAD)
    fhd_arr = np.full(n_cells, np.nan, dtype=np.float32)
    lai_arr = np.full(n_cells, np.nan, dtype=np.float32)
    laisub_arr = np.full(n_cells, np.nan, dtype=np.float32)
    # Reuse the already-sorted groups from GFP/VCI loop
    order2 = np.argsort(cid_clean)
    z_s2 = z_clean[order2]
    cid_s2 = cid_clean[order2]
    if len(cid_s2) > 0:
        changes2 = np.where(np.diff(cid_s2) != 0)[0] + 1
        starts2 = np.concatenate([[0], changes2])
        ends2 = np.concatenate([changes2, [len(z_s2)]])
        gids2 = cid_s2[starts2]
        for i in range(len(starts2)):
            s, e = starts2[i], ends2[i]
            if e - s >= 5:
                fhd_v, lai_v, laisub_v = _fhd_lai_lai_sub(z_s2[s:e])
                fhd_arr[gids2[i]] = fhd_v
                lai_arr[gids2[i]] = lai_v
                laisub_arr[gids2[i]] = laisub_v
    metrics[18] = fhd_arr.reshape(n_rows, n_cols)
    metrics[19] = lai_arr.reshape(n_rows, n_cols)
    metrics[20] = laisub_arr.reshape(n_rows, n_cols)

    # Bands 8-9: mean_sd, sd_sd (3m sub-grid within each 10m cell)
    n_sub = max(1, CELL_SIZE // SUBCELL_SIZE)  # 10/3 → 3 sub-cells per side
    if n_sub > 1:
        sub_total_cols = n_cols * n_sub
        sub_total_rows = n_rows * n_sub

        # Use cleaned x, y, z for subcell assignment
        veg_clean = z_clean >= Z0
        xv = x_clean[veg_clean]
        yv = y_clean[veg_clean]
        zv_sub = z_clean[veg_clean]

        sc = np.clip(((xv - x_min) / SUBCELL_SIZE).astype(np.int32), 0, sub_total_cols - 1)
        sr = np.clip(((y_max - yv) / SUBCELL_SIZE).astype(np.int32), 0, sub_total_rows - 1)
        sub_id = sr * sub_total_cols + sc

        n_sub_cells = sub_total_rows * sub_total_cols
        sub_count = np.bincount(sub_id, minlength=n_sub_cells).astype(np.float64)
        sub_sum = np.bincount(sub_id, weights=zv_sub, minlength=n_sub_cells)
        sub_sumsq = np.bincount(sub_id, weights=zv_sub**2, minlength=n_sub_cells)

        sub_valid = sub_count > 1
        sub_var = np.full(n_sub_cells, np.nan)
        sub_var[sub_valid] = np.maximum(
            (sub_sumsq[sub_valid] - sub_sum[sub_valid]**2 / sub_count[sub_valid]) / (sub_count[sub_valid] - 1), 0)
        sub_sd = np.sqrt(sub_var)

        sub_sd_grid = sub_sd.reshape(sub_total_rows, sub_total_cols)
        sub_valid_grid = sub_valid.reshape(sub_total_rows, sub_total_cols)

        for r in range(n_rows):
            for c in range(n_cols):
                block = sub_sd_grid[r*n_sub:(r+1)*n_sub, c*n_sub:(c+1)*n_sub]
                vmask_b = sub_valid_grid[r*n_sub:(r+1)*n_sub, c*n_sub:(c+1)*n_sub]
                vals = block[vmask_b]
                if len(vals) > 1:
                    metrics[8, r, c] = np.mean(vals)
                    metrics[9, r, c] = np.std(vals, ddof=1)
                elif len(vals) == 1:
                    metrics[8, r, c] = vals[0]
                    metrics[9, r, c] = 0.0

    # CRS
    return metrics, epsg, (x_min, y_min, x_max, y_max)


# ─── Site Processing & Mosaic ────────────────────────────────────────────────

def discover_sites(base_path):
    sites = {}
    for product_dir in ["neon-aop-products", "neon-aop-provisional-products"]:
        product_path = base_path / product_dir
        if not product_path.exists():
            continue
        for year_dir in sorted(product_path.iterdir()):
            if not year_dir.is_dir() or not year_dir.name.isdigit():
                continue
            full_site = year_dir / "FullSite"
            if not full_site.exists():
                continue
            for domain_dir in sorted(full_site.iterdir()):
                if not domain_dir.is_dir():
                    continue
                for site_dir in sorted(domain_dir.iterdir()):
                    if not site_dir.is_dir():
                        continue
                    pc_dir = site_dir / "L1" / "DiscreteLidar" / "ClassifiedPointCloud"
                    if not pc_dir.exists():
                        continue
                    laz_files = sorted(pc_dir.glob("NEON_*.laz"))
                    if laz_files:
                        sites[site_dir.name] = laz_files
    return sites


def process_site(site_key, laz_files, n_workers=1):
    print(f"\n{'='*60}")
    print(f"Processing {site_key}: {len(laz_files)} tiles")
    print(f"{'='*60}")

    t0 = time.time()

    if n_workers > 1:
        results = Parallel(n_jobs=n_workers, backend="loky")(
            delayed(process_tile)(f) for f in tqdm(laz_files, desc=f"  {site_key}")
        )
    else:
        results = [process_tile(f) for f in tqdm(laz_files, desc=f"  {site_key}")]

    valid = [r for r in results if r is not None]
    if not valid:
        print(f"  [WARN] No valid tiles for {site_key}")
        return

    all_bounds = [v[2] for v in valid]
    global_xmin = min(b[0] for b in all_bounds)
    global_ymin = min(b[1] for b in all_bounds)
    global_xmax = max(b[2] for b in all_bounds)
    global_ymax = max(b[3] for b in all_bounds)

    # Snap to Sentinel-2 grid (10m multiples — already snapped by construction)
    global_xmin = np.floor(global_xmin / CELL_SIZE) * CELL_SIZE
    global_ymin = np.floor(global_ymin / CELL_SIZE) * CELL_SIZE
    global_xmax = np.ceil(global_xmax / CELL_SIZE) * CELL_SIZE
    global_ymax = np.ceil(global_ymax / CELL_SIZE) * CELL_SIZE

    total_cols = int((global_xmax - global_xmin) / CELL_SIZE)
    total_rows = int((global_ymax - global_ymin) / CELL_SIZE)

    mosaic = np.full((N_BANDS, total_rows, total_cols), np.nan, dtype=np.float32)

    for grid, epsg_tile, bounds in valid:
        x0, y0, x1, y1 = bounds
        col_off = int((x0 - global_xmin) / CELL_SIZE)
        row_off = int((global_ymax - y1) / CELL_SIZE)
        h, w = grid.shape[1], grid.shape[2]
        r_end = min(row_off + h, total_rows)
        c_end = min(col_off + w, total_cols)
        mosaic[:, row_off:r_end, col_off:c_end] = grid[:, :r_end-row_off, :c_end-col_off]

    # Get EPSG from valid results (use majority vote)
    epsg = None
    epsg_counts = {}
    for _, ep, _ in valid:
        if ep:
            epsg_counts[ep] = epsg_counts.get(ep, 0) + 1
    if epsg_counts:
        epsg = max(epsg_counts, key=epsg_counts.get)
        if len(epsg_counts) > 1:
            print(f"  WARNING: Mixed CRS across tiles: {epsg_counts}. Using EPSG:{epsg}")
    else:
        print(f"  ERROR: No CRS found for any tile in {site_key}. Output will lack CRS!")

    transform = from_origin(global_xmin, global_ymax, CELL_SIZE, CELL_SIZE)
    out_path = OUTPUT_DIR / f"{site_key}_FSD_{CELL_SIZE}m.tif"

    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "width": total_cols,
        "height": total_rows,
        "count": N_BANDS,
        "transform": transform,
        "compress": "lzw",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
        "nodata": np.nan,
    }
    if epsg:
        profile["crs"] = CRS.from_epsg(epsg)

    with rasterio.open(str(out_path), "w", **profile) as dst:
        for i in range(N_BANDS):
            dst.write(mosaic[i], i + 1)
            dst.set_band_description(i + 1, BAND_NAMES[i])

    elapsed = time.time() - t0
    fsize_mb = out_path.stat().st_size / 1024 / 1024
    print(f"  Saved: {out_path} ({fsize_mb:.1f} MB)")
    print(f"  Size: {total_rows} x {total_cols} cells ({total_rows*CELL_SIZE/1000:.1f} x {total_cols*CELL_SIZE/1000:.1f} km)")
    print(f"  CRS: EPSG:{epsg}")
    print(f"  Elapsed: {elapsed:.1f}s ({elapsed/len(laz_files):.1f}s/tile)")


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compute NEON LiDAR structural diversity metrics")
    parser.add_argument("--site", type=str, default=None, help="Filter by site code (e.g., HARV)")
    parser.add_argument("--year", type=str, default=None, help="Filter by year (e.g., 2022)")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers per site")
    parser.add_argument("--dry-run", action="store_true", help="List sites without processing")
    parser.add_argument("--cell-size", type=int, default=None,
                        help="Override cell size in meters (default 10). "
                             "Results saved to structural_diversity_{N}m/")
    args = parser.parse_args()

    # Allow overriding cell size (e.g., 1m for plot-level analysis)
    global CELL_SIZE, SUBCELL_SIZE, OUTPUT_DIR
    if args.cell_size is not None:
        CELL_SIZE = args.cell_size
        SUBCELL_SIZE = max(1, CELL_SIZE // 3) if CELL_SIZE >= 3 else 1
        OUTPUT_DIR = Path(f"E:/neon_lidar/structural_diversity_{CELL_SIZE}m")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Discovering NEON sites...")
    sites = discover_sites(NEON_BASE)
    print(f"Found {len(sites)} site-year combinations, {sum(len(v) for v in sites.values())} total tiles")

    if args.site:
        sites = {k: v for k, v in sites.items() if args.site.upper() in k.upper()}
    if args.year:
        sites = {k: v for k, v in sites.items() if k.startswith(args.year)}

    if not sites:
        print("No matching sites found.")
        return

    suffix = f"_FSD_{CELL_SIZE}m"
    existing = {f.stem.replace(suffix, "") for f in OUTPUT_DIR.glob(f"*{suffix}.tif")}
    to_process = {k: v for k, v in sites.items() if k not in existing}
    skipped = len(sites) - len(to_process)
    if skipped > 0:
        print(f"Skipping {skipped} already-processed sites")

    print(f"\nWill process {len(to_process)} site-years:")
    for k in sorted(to_process.keys()):
        print(f"  {k}: {len(to_process[k])} tiles")

    if args.dry_run:
        return

    for site_key in sorted(to_process.keys()):
        process_site(site_key, to_process[site_key], n_workers=args.workers)

    print(f"\n{'='*60}")
    print(f"All done! Results in: {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
