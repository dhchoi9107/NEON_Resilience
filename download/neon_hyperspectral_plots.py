"""
NEON Hyperspectral Download — Plot Buffers Only
=================================================
Downloads NEON AOP hyperspectral reflectance (DP3.30006.001) tiles
only for 1km tiles that overlap with vegetation plot locations (± 50m buffer).

This avoids downloading full-site mosaics (10-50x less data).

Product: DP3.30006.001 — Spectrometer orthorectified surface directional
         reflectance — mosaic (426 bands, 1m resolution, HDF5 format)

Usage:
  python download/neon_hyperspectral_plots.py estimate       # check availability & size
  python download/neon_hyperspectral_plots.py download        # download all
  python download/neon_hyperspectral_plots.py download --site HARV  # single site
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import os
import re
import json
import time
import argparse
import urllib.request
import urllib.error

from site_config import SITES, YEARS, VEG_STRUCT_DIR, TOKEN, NEON_API, api_get

import pandas as pd

DPID = "DP3.30006.002"  # bidirectional (BRDF) corrected mosaic
BUFFER = 50  # meters around plot edge
OUTPUT_DIR = Path("E:/neon_lidar/hyperspectral_plots")


def load_plot_tiles():
    """Determine which 1km tiles are needed per site based on plot locations.

    Returns:
      site_tiles: {site: set of (easting, northing)}
      plot_coords: DataFrame with siteID, plotID, easting, northing
    """
    ppy = pd.read_csv(VEG_STRUCT_DIR / "vst_perplotperyear.csv", low_memory=False)
    coords = ppy[["siteID", "plotID", "easting", "northing"]].dropna()
    coords = coords.drop_duplicates(subset=["siteID", "plotID"])

    site_tiles = {}  # {site: set of (easting, northing)}
    for _, row in coords.iterrows():
        site = row["siteID"]
        px, py = row["easting"], row["northing"]
        if site not in site_tiles:
            site_tiles[site] = set()
        for te in range(int((px - BUFFER) // 1000) * 1000,
                        int((px + BUFFER) // 1000) * 1000 + 1001, 1000):
            for tn in range(int((py - BUFFER) // 1000) * 1000,
                            int((py + BUFFER) // 1000) * 1000 + 1001, 1000):
                site_tiles[site].add((te, tn))

    return site_tiles, coords


def get_available_files(site, month):
    """Get file list from NEON API for a site-month."""
    url = f"{NEON_API}/data/{DPID}/{site}/{month}"
    try:
        data = api_get(url)
        return data.get("data", {}).get("files", [])
    except Exception:
        return []


def tile_matches(filename, needed_tiles):
    """Check if a filename corresponds to a needed tile."""
    m = re.search(r'_(\d{6})_(\d{7})_', filename)
    if not m:
        return False
    te, tn = int(m.group(1)), int(m.group(2))
    return (te, tn) in needed_tiles


def download_file(url, dest_path, token=None):
    """Download a file with resume support."""
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    if dest_path.exists():
        return "skip"

    req = urllib.request.Request(url)
    if token:
        req.add_header("X-API-Token", token)

    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                with open(str(dest_path), "wb") as f:
                    while True:
                        chunk = r.read(1024 * 1024)  # 1MB chunks
                        if not chunk:
                            break
                        f.write(chunk)
            return "ok"
        except Exception as e:
            if attempt < 2:
                time.sleep(5)
            else:
                return f"error: {e}"


def estimate(site_filter=None, latest_only=False):
    """Check availability and estimate download size."""
    site_tiles, _ = load_plot_tiles()
    sites = [site_filter] if site_filter else [s for s in SITES if s in site_tiles]

    print(f"Product: {DPID} (Hyperspectral Reflectance - Mosaic, BRDF corrected)")
    print(f"Strategy: Only tiles overlapping plot locations (± {BUFFER}m buffer)")
    print(f"Output: {OUTPUT_DIR}\n")

    # Get availability
    prod_data = api_get(f"{NEON_API}/products/{DPID}")["data"]
    avail = {}
    for s in prod_data["siteCodes"]:
        avail[s["siteCode"]] = sorted(s.get("availableMonths", []))

    total_tiles = 0
    total_files = 0
    total_size = 0

    for site in sites:
        months = avail.get(site, [])
        year_months = [m for m in months if any(m.startswith(str(y)) for y in YEARS)]
        if latest_only and year_months:
            latest_yr = year_months[-1][:4]
            year_months = [m for m in year_months if m.startswith(latest_yr)]
        needed = site_tiles.get(site, set())

        if not year_months or not needed:
            print(f"  {site}: no data or no plots")
            continue

        # Check most recent month to estimate file sizes
        recent = year_months[-1]
        files = get_available_files(site, recent)
        matching = [f for f in files
                    if f["name"].endswith(".h5") and tile_matches(f["name"], needed)]

        n_tiles = len(needed)
        n_files = len(matching)
        size_mb = sum(f.get("size", 0) for f in matching) / 1e6

        # Extrapolate to all available years
        n_years = len(set(m[:4] for m in year_months))
        est_total_mb = size_mb * n_years

        total_tiles += n_tiles
        total_files += n_files * n_years
        total_size += est_total_mb

        print(f"  {site}: {n_tiles:3d} tiles, {n_years} years, "
              f"~{n_files} files/year, ~{size_mb:.0f} MB/year, "
              f"est total: {est_total_mb/1000:.1f} GB")

    print(f"\n  TOTAL: ~{total_tiles} unique tiles, ~{total_files} files, "
          f"~{total_size/1000:.0f} GB estimated")


def clip_h5_to_plots(h5_path, site_plots, site_dir, year_str="unknown"):
    """Download full 1km tile, clip to each overlapping plot (100x100m), save as GeoTIFF, delete H5.

    For each plot overlapping this tile, extracts a 100x100m window (plot center +/- 50m)
    and saves as a multi-band GeoTIFF (~426 bands). Much smaller than full H5.
    """
    import h5py
    import numpy as np

    m = re.search(r'_(\d{6})_(\d{7})_', h5_path.name)
    if not m:
        return 0
    tile_e, tile_n = int(m.group(1)), int(m.group(2))

    # Find plots that overlap with this tile (within BUFFER)
    overlapping = []
    for _, row in site_plots.iterrows():
        px, py = row["easting"], row["northing"]
        if (tile_e <= px + BUFFER and tile_e + 1000 >= px - BUFFER and
                tile_n <= py + BUFFER and tile_n + 1000 >= py - BUFFER):
            overlapping.append(row)

    if not overlapping:
        h5_path.unlink(missing_ok=True)
        return 0

    n_clipped = 0
    try:
        with h5py.File(str(h5_path), "r") as hf:
            # Navigate HDF5 structure to find reflectance data
            # Typical: /{site}/Reflectance/Reflectance_Data
            site_key = list(hf.keys())[0]
            refl_grp = hf[site_key]["Reflectance"]
            refl_data = refl_grp["Reflectance_Data"]  # shape: (bands, rows, cols)

            # Get spatial info
            map_info = refl_grp["Metadata"]["Coordinate_System"]["Map_Info"][()].decode()
            parts = map_info.split(",")
            # Map_Info: UTM, 1, 1, easting, northing, xres, yres, zone, ...
            x_origin = float(parts[3])
            y_origin = float(parts[4])
            x_res = float(parts[5])
            y_res = float(parts[6])

            n_bands = refl_data.shape[0]

            for prow in overlapping:
                px, py = prow["easting"], prow["northing"]
                plot_id = prow["plotID"]

                # Window: plot center +/- 50m
                xmin = px - BUFFER
                xmax = px + BUFFER
                ymin = py - BUFFER
                ymax = py + BUFFER

                col_start = max(0, int((xmin - x_origin) / x_res))
                col_end = min(refl_data.shape[2], int((xmax - x_origin) / x_res))
                row_start = max(0, int((y_origin - ymax) / y_res))
                row_end = min(refl_data.shape[1], int((y_origin - ymin) / y_res))

                if col_end <= col_start or row_end <= row_start:
                    continue

                # Read only the window (memory efficient)
                window_data = refl_data[:, row_start:row_end, col_start:col_end]

                # Save as multi-band GeoTIFF
                import rasterio
                from rasterio.crs import CRS
                from rasterio.transform import from_origin
                from site_config import SITE_EPSG

                out_name = f"{year_str}_{plot_id}_hyper.tif"
                out_path = site_dir / out_name
                if out_path.exists():
                    n_clipped += 1
                    continue

                transform = from_origin(
                    x_origin + col_start * x_res,
                    y_origin - row_start * y_res,
                    x_res, y_res
                )
                epsg = SITE_EPSG.get(prow["siteID"])
                profile = {
                    "driver": "GTiff",
                    "dtype": "int16",
                    "width": col_end - col_start,
                    "height": row_end - row_start,
                    "count": n_bands,
                    "transform": transform,
                    "compress": "lzw",
                    "tiled": True,
                    "nodata": -9999,
                }
                if epsg:
                    profile["crs"] = CRS.from_epsg(epsg)

                with rasterio.open(str(out_path), "w", **profile) as dst:
                    for b in range(n_bands):
                        dst.write(window_data[b].astype(np.int16), b + 1)

                n_clipped += 1

    except Exception as e:
        print(f"      Clip error: {e}")

    # Delete original H5
    h5_path.unlink(missing_ok=True)
    return n_clipped


def download(site_filter=None, latest_only=False):
    """Download hyperspectral tiles, clip to plot windows, delete originals."""
    site_tiles, plot_coords = load_plot_tiles()
    sites = [site_filter] if site_filter else [s for s in SITES if s in site_tiles]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Get availability
    prod_data = api_get(f"{NEON_API}/products/{DPID}")["data"]
    avail = {}
    for s in prod_data["siteCodes"]:
        avail[s["siteCode"]] = sorted(s.get("availableMonths", []))

    for site in sites:
        months = avail.get(site, [])
        year_months = [m for m in months if any(m.startswith(str(y)) for y in YEARS)]
        if latest_only and year_months:
            latest_yr = year_months[-1][:4]
            year_months = [m for m in year_months if m.startswith(latest_yr)]
        needed = site_tiles.get(site, set())
        site_plots = plot_coords[plot_coords["siteID"] == site]

        if not year_months or not needed:
            continue

        site_dir = OUTPUT_DIR / site
        site_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"  {site}: {len(needed)} tiles, {len(year_months)} months, {len(site_plots)} plots")

        for month in year_months:
            marker = site_dir / f".done_{month}"
            if marker.exists():
                continue

            files = get_available_files(site, month)
            matching = [f for f in files
                        if f["name"].endswith(".h5") and tile_matches(f["name"], needed)]

            if not matching:
                continue

            # Determine which tiles actually need downloading:
            # skip if ALL plots for that tile are already clipped
            month_year = month[:4]
            tiles_to_download = []
            for finfo in matching:
                tm = re.search(r'_(\d{6})_(\d{7})_', finfo["name"])
                if not tm:
                    continue
                te, tn = int(tm.group(1)), int(tm.group(2))

                # Which plots overlap this tile?
                plots_for_tile = []
                for _, prow in site_plots.iterrows():
                    px, py = prow["easting"], prow["northing"]
                    if (te <= px + BUFFER and te + 1000 >= px - BUFFER and
                            tn <= py + BUFFER and tn + 1000 >= py - BUFFER):
                        plots_for_tile.append(prow["plotID"])

                # Check if all clips already exist
                all_exist = all(
                    (site_dir / f"{month_year}_{pid}_hyper.tif").exists()
                    for pid in plots_for_tile
                )
                if all_exist and plots_for_tile:
                    continue  # skip this tile entirely

                tiles_to_download.append(finfo)

            if not tiles_to_download:
                marker.touch()
                continue

            n_clipped_total = 0
            for finfo in tiles_to_download:
                tmp_path = site_dir / f"_tmp_{finfo['name']}"

                result = download_file(finfo["url"], tmp_path, TOKEN)
                if result == "ok":
                    n = clip_h5_to_plots(tmp_path, site_plots, site_dir, year_str=month[:4])
                    n_clipped_total += n
                elif result == "skip":
                    if tmp_path.exists():
                        n = clip_h5_to_plots(tmp_path, site_plots, site_dir, year_str=month[:4])
                        n_clipped_total += n
                else:
                    print(f"    FAIL: {finfo['name']}: {result}")
                    tmp_path.unlink(missing_ok=True)

            print(f"    {month}: {n_clipped_total} plot clips from {len(matching)} tiles")

            marker.touch()


def main():
    parser = argparse.ArgumentParser(
        description="Download NEON hyperspectral tiles for plot locations only")
    parser.add_argument("mode", choices=["estimate", "download"],
                        help="estimate: check size; download: fetch data")
    parser.add_argument("--site", type=str, default=None,
                        help="Process single site (e.g., HARV)")
    parser.add_argument("--latest-only", action="store_true",
                        help="Download only the most recent year per site")
    args = parser.parse_args()

    if args.mode == "estimate":
        estimate(site_filter=args.site, latest_only=args.latest_only)
    else:
        download(site_filter=args.site, latest_only=args.latest_only)


if __name__ == "__main__":
    main()
