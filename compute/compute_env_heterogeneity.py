"""
Environmental Heterogeneity Metrics
=====================================
Computes three classes of environmental heterogeneity per site-year:

1) Compositional heterogeneity
   - CV and range of canopy height (from FSD rasters)
   - CV of NDVI (from spectral diversity rasters, if available)

2) Configurational heterogeneity (gradient surface metrics)
   - Sa (average roughness), Sq (RMS roughness)
   - Ssk (surface skewness), Sku (surface kurtosis)
   - Computed at multiple window sizes (100m, 500m, 1km)

3) Fragmentation (landscape metrics from NLCD)
   - Shannon Landscape Diversity Index (SHDI)
   - Patch Density (PD)
   - Number of land cover classes

Usage:
  python compute_env_heterogeneity.py                  # all sites
  python compute_env_heterogeneity.py --site HARV      # single site
  python compute_env_heterogeneity.py --skip-nlcd      # skip NLCD download
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import warnings

import numpy as np
import rasterio
from scipy.ndimage import uniform_filter, label

from site_config import (
    SITES, SITE_EPSG, SITE_COORDS, FSD_DIR, SPEC_DIV_DIR, ENV_HET_DIR,
    CELL_SIZE, get_fsd_files, get_fsd_band_index,
)

warnings.filterwarnings("ignore")

CHM_BAND = get_fsd_band_index("mean_max_canopy_ht")
WINDOW_SIZES_M = [100, 500, 1000]  # in meters


# ─── 1) Compositional Heterogeneity ─────────────────────────────────────────

def compositional_heterogeneity(fsd_path):
    """Compute site-level compositional heterogeneity from canopy height."""
    with rasterio.open(str(fsd_path)) as src:
        chm = src.read(CHM_BAND).astype(np.float32)
        chm[chm == src.nodata] = np.nan if src.nodata else chm[chm < 0]

    valid = chm[np.isfinite(chm)]
    if len(valid) < 10:
        return {}

    mean_h = float(np.mean(valid))
    sd_h = float(np.std(valid))
    cv_h = sd_h / mean_h if mean_h > 0 else np.nan

    return {
        "chm_mean": mean_h,
        "chm_sd": sd_h,
        "chm_cv": cv_h,
        "chm_range": float(np.ptp(valid)),
        "chm_iqr": float(np.percentile(valid, 75) - np.percentile(valid, 25)),
    }


# ─── 2) Configurational Heterogeneity (Gradient Surface Metrics) ────────────

def gradient_surface_metrics(fsd_path, window_sizes_m=WINDOW_SIZES_M):
    """Compute gradient surface metrics at multiple window sizes.

    Metrics (computed per window, then site-level mean):
      Sa  = mean(|z - z_mean_local|)       (average roughness)
      Sq  = sqrt(mean((z - z_mean_local)^2)) (RMS roughness)
      Ssk = mean((z - z_mean)^3) / Sq^3    (skewness)
      Sku = mean((z - z_mean)^4) / Sq^4    (kurtosis)
    """
    with rasterio.open(str(fsd_path)) as src:
        chm = src.read(CHM_BAND).astype(np.float64)

    chm[~np.isfinite(chm)] = np.nan

    results = {}
    for win_m in window_sizes_m:
        win_px = max(win_m // CELL_SIZE, 3)  # window in pixels
        if win_px % 2 == 0:
            win_px += 1  # ensure odd

        # Local mean via uniform filter (handles NaN by replacing temporarily)
        filled = np.where(np.isfinite(chm), chm, 0.0)
        count = np.where(np.isfinite(chm), 1.0, 0.0)

        sum_local = uniform_filter(filled, size=win_px, mode="constant", cval=0.0)
        cnt_local = uniform_filter(count, size=win_px, mode="constant", cval=0.0)
        cnt_local[cnt_local < 1] = np.nan
        mean_local = sum_local / cnt_local

        deviation = chm - mean_local
        dev_valid = deviation[np.isfinite(deviation)]

        if len(dev_valid) < 10:
            continue

        sa = float(np.mean(np.abs(dev_valid)))
        sq = float(np.sqrt(np.mean(dev_valid ** 2)))

        ssk = np.nan
        sku = np.nan
        if sq > 1e-6:
            ssk = float(np.mean(dev_valid ** 3) / sq ** 3)
            sku = float(np.mean(dev_valid ** 4) / sq ** 4)

        suffix = f"_{win_m}m"
        results[f"Sa{suffix}"] = sa
        results[f"Sq{suffix}"] = sq
        results[f"Ssk{suffix}"] = ssk
        results[f"Sku{suffix}"] = sku

    return results


# ─── 3) Fragmentation (NLCD) ────────────────────────────────────────────────

def download_nlcd_clip(site, buffer_m=5000, nlcd_year=2021):
    """Download NLCD clip for a site from MRLC WCS."""
    import urllib.request
    import urllib.error

    if site not in SITE_COORDS or site not in SITE_EPSG:
        return None

    lat, lon = SITE_COORDS[site]
    epsg = SITE_EPSG[site]

    # Approximate UTM center from site coordinates
    # Use the FSD raster bounds instead for accurate UTM extent
    fsd_entries = get_fsd_files(site)
    if not fsd_entries:
        return None

    # Use the most recent FSD file for bounds
    _, _, fsd_path = fsd_entries[-1]
    with rasterio.open(str(fsd_path)) as src:
        bounds = src.bounds
        crs = src.crs

    xmin = bounds.left - buffer_m
    ymin = bounds.bottom - buffer_m
    xmax = bounds.right + buffer_m
    ymax = bounds.top + buffer_m

    # MRLC WCS request
    width = int((xmax - xmin) / 30)
    height = int((ymax - ymin) / 30)

    wcs_url = (
        f"https://www.mrlc.gov/geoserver/mrlc_display/"
        f"NLCD_{nlcd_year}_Land_Cover_L48/wcs?"
        f"service=WCS&version=2.0.1&request=GetCoverage"
        f"&CoverageId=NLCD_{nlcd_year}_Land_Cover_L48"
        f"&subset=X({xmin},{xmax})"
        f"&subset=Y({ymin},{ymax})"
        f"&subsettingCrs=http://www.opengis.net/def/crs/EPSG/0/{epsg}"
        f"&format=image/tiff"
    )

    out_path = ENV_HET_DIR / f"{site}_NLCD_{nlcd_year}.tif"
    if out_path.exists():
        return out_path

    print(f"    Downloading NLCD for {site}...")
    try:
        urllib.request.urlretrieve(wcs_url, str(out_path))
        # Verify it's a valid GeoTIFF
        with rasterio.open(str(out_path)) as src:
            if src.width < 10 or src.height < 10:
                out_path.unlink()
                return None
        print(f"    Saved: {out_path.name} ({width}x{height})")
        return out_path
    except Exception as e:
        print(f"    NLCD download failed for {site}: {e}")
        if out_path.exists():
            out_path.unlink()
        return None


def fragmentation_metrics(nlcd_path):
    """Compute landscape fragmentation metrics from NLCD raster."""
    with rasterio.open(str(nlcd_path)) as src:
        lc = src.read(1)

    # NLCD classes: 11=water, 21-24=developed, 31=barren, 41-43=forest,
    # 52=shrub, 71=grass, 81-82=agriculture, 90-95=wetland
    # Exclude nodata (0 or 255)
    valid = lc[(lc > 0) & (lc < 255)]
    if len(valid) < 100:
        return {}

    classes, counts = np.unique(valid, return_counts=True)
    n_classes = len(classes)
    total = counts.sum()

    # Shannon Landscape Diversity Index
    proportions = counts / total
    shdi = -np.sum(proportions * np.log(proportions))

    # Patch density: count distinct patches
    total_patches = 0
    for cls in classes:
        mask = (lc == cls).astype(int)
        labeled, n_patches = label(mask)
        total_patches += n_patches

    # Area in km^2 (30m pixels)
    area_km2 = total * 30 * 30 / 1e6
    patch_density = total_patches / area_km2 if area_km2 > 0 else np.nan

    # Forest proportion
    forest_mask = np.isin(valid, [41, 42, 43, 90])
    forest_prop = float(np.sum(forest_mask) / len(valid))

    return {
        "nlcd_n_classes": n_classes,
        "nlcd_shdi": float(shdi),
        "nlcd_patch_density": float(patch_density),
        "nlcd_forest_proportion": forest_prop,
        "nlcd_area_km2": float(area_km2),
    }


# ─── Main Processing ────────────────────────────────────────────────────────

def process_site(site, skip_nlcd=False):
    """Compute all heterogeneity metrics for a site."""
    fsd_entries = get_fsd_files(site)
    if not fsd_entries:
        print(f"  SKIP {site}: no FSD files")
        return []

    print(f"  {site}: {len(fsd_entries)} site-years")

    results = []
    for year, st, fsd_path in fsd_entries:
        row = {"siteID": site, "year": year}

        # 1) Compositional
        comp = compositional_heterogeneity(fsd_path)
        row.update(comp)

        # 2) Configurational (gradient surface metrics)
        grad = gradient_surface_metrics(fsd_path)
        row.update(grad)

        results.append(row)

    # 3) Fragmentation (NLCD) - once per site (static land cover)
    nlcd_metrics = {}
    if not skip_nlcd:
        nlcd_path = download_nlcd_clip(site)
        if nlcd_path:
            nlcd_metrics = fragmentation_metrics(nlcd_path)

    # Append NLCD metrics to all years for this site
    for row in results:
        row.update(nlcd_metrics)

    if results:
        print(f"    CHM CV: {results[-1].get('chm_cv', 'N/A'):.3f}, "
              f"Sa_500m: {results[-1].get('Sa_500m', 'N/A'):.2f}")
        if nlcd_metrics:
            print(f"    NLCD SHDI: {nlcd_metrics.get('nlcd_shdi', 'N/A'):.3f}, "
                  f"Forest: {nlcd_metrics.get('nlcd_forest_proportion', 0):.1%}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Compute environmental heterogeneity")
    parser.add_argument("--site", type=str, default=None, help="Process single site")
    parser.add_argument("--skip-nlcd", action="store_true", help="Skip NLCD download")
    args = parser.parse_args()

    ENV_HET_DIR.mkdir(parents=True, exist_ok=True)

    sites = [args.site] if args.site else SITES
    print(f"Environmental Heterogeneity Computation")
    print(f"Output: {ENV_HET_DIR}")
    print(f"Sites: {len(sites)}\n")

    all_results = []
    for site in sites:
        print(f"\n{'='*50}")
        rows = process_site(site, skip_nlcd=args.skip_nlcd)
        all_results.extend(rows)

    if all_results:
        import pandas as pd
        df = pd.DataFrame(all_results)
        csv_path = ENV_HET_DIR / "heterogeneity_all.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path} ({len(df)} site-years)")


if __name__ == "__main__":
    main()
