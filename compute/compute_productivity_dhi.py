"""
Dynamic Habitat Indices (DHI) from Multi-temporal Canopy Height
================================================================
Computes productivity proxies from inter-annual canopy height changes
derived from existing FSD (Forest Structural Diversity) rasters.

Metrics (per pixel, 10m):
  Cumulative DHI  - sum of positive height changes (total growth)
  Minimum DHI     - minimum annual height change (environmental stress)
  Variation DHI   - SD of annual height changes (stability)
  Height trend    - linear slope of height vs year (long-term growth rate)

Usage:
  python compute_productivity_dhi.py                    # all sites
  python compute_productivity_dhi.py --site HARV        # single site
  python compute_productivity_dhi.py --dry-run          # list only
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import warnings
from collections import defaultdict

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin

from site_config import (
    SITES, SITE_EPSG, FSD_DIR, DHI_DIR, CELL_SIZE,
    get_fsd_files, get_fsd_band_index, FSD_BANDS,
)

warnings.filterwarnings("ignore")

# Band index for canopy height metric (1-based)
CHM_BAND = get_fsd_band_index("mean_max_canopy_ht")


def compute_site_dhi(site):
    """Compute DHI metrics for a single site across all available years."""

    # Collect FSD files for this site
    fsd_entries = get_fsd_files(site)
    if len(fsd_entries) < 2:
        print(f"  SKIP {site}: need >= 2 years, found {len(fsd_entries)}")
        return None

    years = [yr for yr, _, _ in fsd_entries]
    paths = [p for _, _, p in fsd_entries]
    print(f"  {site}: {len(years)} years ({min(years)}-{max(years)})")

    # ── Read metadata from all files, find common extent ──
    profiles = []
    for p in paths:
        with rasterio.open(str(p)) as src:
            profiles.append({
                "bounds": src.bounds,
                "transform": src.transform,
                "width": src.width,
                "height": src.height,
                "crs": src.crs,
            })

    # Intersection bounding box (common area across all years)
    xmin = max(p["bounds"].left for p in profiles)
    ymin = max(p["bounds"].bottom for p in profiles)
    xmax = min(p["bounds"].right for p in profiles)
    ymax = min(p["bounds"].top for p in profiles)

    if xmax <= xmin or ymax <= ymin:
        print(f"  SKIP {site}: no spatial overlap across years")
        return None

    # Snap to grid
    xmin = np.floor(xmin / CELL_SIZE) * CELL_SIZE
    ymin = np.floor(ymin / CELL_SIZE) * CELL_SIZE
    xmax = np.ceil(xmax / CELL_SIZE) * CELL_SIZE
    ymax = np.ceil(ymax / CELL_SIZE) * CELL_SIZE

    ncols = int((xmax - xmin) / CELL_SIZE)
    nrows = int((ymax - ymin) / CELL_SIZE)

    if ncols < 1 or nrows < 1:
        print(f"  SKIP {site}: overlap too small ({ncols}x{nrows})")
        return None

    # ── Read CHM band from each year into a 3D stack ──
    stack = np.full((len(years), nrows, ncols), np.nan, dtype=np.float32)

    for i, p in enumerate(paths):
        with rasterio.open(str(p)) as src:
            window = rasterio.windows.from_bounds(
                xmin, ymin, xmax, ymax, transform=src.transform
            )
            # Clip window to source bounds
            window = window.intersection(
                rasterio.windows.Window(0, 0, src.width, src.height)
            )
            data = src.read(CHM_BAND, window=window)

            # Place into stack (handle minor size mismatches from rounding)
            h, w = data.shape
            h = min(h, nrows)
            w = min(w, ncols)
            stack[i, :h, :w] = data[:h, :w]

    # ── Compute annual height changes ──
    # For non-consecutive years, normalize to annual rate
    year_arr = np.array(years, dtype=np.float32)
    n_intervals = len(years) - 1

    if n_intervals < 1:
        print(f"  SKIP {site}: need at least 2 years for DHI")
        return None

    # Annual height change rates: delta_h / delta_year
    delta_h = np.full((n_intervals, nrows, ncols), np.nan, dtype=np.float32)
    for i in range(n_intervals):
        dt = year_arr[i + 1] - year_arr[i]
        if dt > 0:
            delta_h[i] = (stack[i + 1] - stack[i]) / dt

    # ── DHI metrics ──
    # Mask: need at least 2 valid intervals for meaningful statistics
    valid_count = np.sum(np.isfinite(delta_h), axis=0)

    # Cumulative DHI: sum of positive height changes (total growth)
    pos_delta = np.where(delta_h > 0, delta_h, 0)
    cumulative = np.nansum(pos_delta, axis=0)
    cumulative[valid_count < 1] = np.nan

    # Minimum DHI: minimum annual change (stress indicator)
    with np.errstate(all="ignore"):
        minimum = np.nanmin(delta_h, axis=0)
    minimum[valid_count < 1] = np.nan

    # Variation DHI: SD of annual changes (stability)
    with np.errstate(all="ignore"):
        variation = np.nanstd(delta_h, axis=0, ddof=0)
    variation[valid_count < 2] = np.nan

    # ── Pixel-wise linear trend (height ~ year) ──
    # slope = cov(year, height) / var(year)
    valid_mask = np.isfinite(stack)
    n_valid = np.sum(valid_mask, axis=0).astype(np.float32)

    trend = np.full((nrows, ncols), np.nan, dtype=np.float32)
    mask_enough = n_valid >= 3  # need >= 3 points for meaningful trend

    if np.any(mask_enough):
        # Vectorized linear regression
        yr_broadcast = np.broadcast_to(
            year_arr[:, np.newaxis, np.newaxis], stack.shape
        ).copy()
        yr_broadcast[~valid_mask] = np.nan

        with np.errstate(all="ignore"):
            yr_mean = np.nanmean(yr_broadcast, axis=0)
            ht_mean = np.nanmean(stack, axis=0)

            yr_dev = yr_broadcast - yr_mean[np.newaxis, :, :]
            ht_dev = stack - ht_mean[np.newaxis, :, :]
            yr_dev[~valid_mask] = 0
            ht_dev[~valid_mask] = 0

            cov_xy = np.nansum(yr_dev * ht_dev, axis=0)
            var_x = np.nansum(yr_dev ** 2, axis=0)

            slope = np.where(var_x > 0, cov_xy / var_x, np.nan)
            trend[mask_enough] = slope[mask_enough]

    # ── Write outputs ──
    DHI_DIR.mkdir(parents=True, exist_ok=True)
    transform = from_origin(xmin, ymax, CELL_SIZE, CELL_SIZE)
    epsg = SITE_EPSG.get(site)

    profile_base = {
        "driver": "GTiff",
        "dtype": "float32",
        "width": ncols,
        "height": nrows,
        "transform": transform,
        "compress": "lzw",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
        "nodata": np.nan,
    }
    if epsg:
        profile_base["crs"] = CRS.from_epsg(epsg)

    # DHI raster (3 bands)
    dhi_path = DHI_DIR / f"{site}_dhi_10m.tif"
    profile = {**profile_base, "count": 3}
    with rasterio.open(str(dhi_path), "w", **profile) as dst:
        dst.write(cumulative, 1)
        dst.set_band_description(1, "cumulative_dhi")
        dst.write(minimum, 2)
        dst.set_band_description(2, "minimum_dhi")
        dst.write(variation, 3)
        dst.set_band_description(3, "variation_dhi")

    # Trend raster (1 band)
    trend_path = DHI_DIR / f"{site}_trend_10m.tif"
    profile = {**profile_base, "count": 1}
    with rasterio.open(str(trend_path), "w", **profile) as dst:
        dst.write(trend, 1)
        dst.set_band_description(1, "height_trend_m_per_yr")

    # Summary statistics
    summary = {
        "site": site,
        "n_years": len(years),
        "year_min": min(years),
        "year_max": max(years),
        "n_pixels": int(np.sum(np.isfinite(cumulative))),
        "cumulative_mean": float(np.nanmean(cumulative)),
        "cumulative_sd": float(np.nanstd(cumulative)),
        "minimum_mean": float(np.nanmean(minimum)),
        "variation_mean": float(np.nanmean(variation)),
        "trend_mean": float(np.nanmean(trend)),
        "trend_sd": float(np.nanstd(trend)),
    }

    print(f"    DHI saved: {dhi_path.name}  ({nrows}x{ncols}, {len(years)} years)")
    print(f"    Trend: {summary['trend_mean']:.4f} m/yr (mean)")
    print(f"    Cumulative growth: {summary['cumulative_mean']:.2f} m (mean)")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Compute DHI from multi-temporal CHM")
    parser.add_argument("--site", type=str, default=None, help="Process single site")
    parser.add_argument("--dry-run", action="store_true", help="List sites only")
    args = parser.parse_args()

    sites = [args.site] if args.site else SITES

    if args.dry_run:
        for site in sites:
            entries = get_fsd_files(site)
            years = [yr for yr, _, _ in entries]
            print(f"  {site}: {len(years)} years "
                  f"({min(years) if years else 'N/A'}-{max(years) if years else 'N/A'})")
        return

    print(f"Computing DHI for {len(sites)} sites")
    print(f"Output: {DHI_DIR}\n")

    summaries = []
    for site in sites:
        print(f"\n{'='*50}")
        result = compute_site_dhi(site)
        if result:
            summaries.append(result)

    # Write summary CSV
    if summaries:
        import csv
        csv_path = DHI_DIR / "productivity_summary.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
            writer.writeheader()
            writer.writerows(summaries)
        print(f"\nSummary: {csv_path} ({len(summaries)} sites)")


if __name__ == "__main__":
    main()
