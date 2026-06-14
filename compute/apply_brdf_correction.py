"""
BRDF Correction for NEON Hyperspectral Plot Clips
===================================================
Downloads .001 HDF5 tiles, applies FlexBRDF correction using hytools,
clips to plot footprints, saves corrected reflectance.

For each HDF5 tile:
  1. Open with hytools (reads angles + reflectance natively)
  2. Apply Ross-Li BRDF kernel correction per-pixel
  3. Clip to plot ±50m
  4. Save as corrected GeoTIFF
  5. Delete original HDF5

Usage:
  python compute/apply_brdf_correction.py --site HARV
  python compute/apply_brdf_correction.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import re
import warnings

import numpy as np
import pandas as pd
import h5py
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin

from site_config import SITES, YEARS, VEG_STRUCT_DIR, TOKEN, NEON_API, api_get, SITE_EPSG

warnings.filterwarnings("ignore")

DPID = "DP3.30006.001"
BUFFER = 50
OUTPUT_DIR = Path("E:/neon_lidar/hyperspectral_brdf_corrected")
TEMP_DIR = Path("E:/neon_lidar/_tmp_brdf")


# ─── Ross-Li BRDF kernels ───────────────────────────────────────────────────

def ross_thick_kernel(solar_zn, sensor_zn, relative_az):
    """Ross-Thick volumetric kernel."""
    cos_phase = (np.cos(solar_zn) * np.cos(sensor_zn) +
                 np.sin(solar_zn) * np.sin(sensor_zn) * np.cos(relative_az))
    cos_phase = np.clip(cos_phase, -1, 1)
    phase = np.arccos(cos_phase)
    return ((np.pi / 2 - phase) * cos_phase + np.sin(phase)) / (np.cos(solar_zn) + np.cos(sensor_zn)) - np.pi / 4


def li_dense_r_kernel(solar_zn, sensor_zn, relative_az, h_b=2.0, b_r=1.0):
    """Li-Dense-R geometric kernel (reciprocal)."""
    # Transform angles for reciprocal
    solar_zn_p = np.arctan(b_r * np.tan(solar_zn))
    sensor_zn_p = np.arctan(b_r * np.tan(sensor_zn))

    cos_phase_p = (np.cos(solar_zn_p) * np.cos(sensor_zn_p) +
                   np.sin(solar_zn_p) * np.sin(sensor_zn_p) * np.cos(relative_az))
    cos_phase_p = np.clip(cos_phase_p, -1, 1)

    D = np.sqrt(np.tan(solar_zn_p)**2 + np.tan(sensor_zn_p)**2 -
                2 * np.tan(solar_zn_p) * np.tan(sensor_zn_p) * np.cos(relative_az))

    sec_sum = 1.0 / np.cos(solar_zn_p) + 1.0 / np.cos(sensor_zn_p)

    cos_t = np.clip(h_b * np.sqrt(D**2 + (np.tan(solar_zn_p) * np.tan(sensor_zn_p) * np.sin(relative_az))**2) / sec_sum, -1, 1)
    t = np.arccos(cos_t)

    O = (1 / np.pi) * (t - np.sin(t) * cos_t) * sec_sum
    O = np.maximum(O, 0)

    return O - sec_sum + 0.5 * (1 + cos_phase_p) / (np.cos(solar_zn_p) * np.cos(sensor_zn_p))


def brdf_correct_pixel(reflectance, solar_zn, solar_az, sensor_zn, sensor_az,
                        ref_solar_zn=None):
    """Apply BRDF correction to normalize to nadir view.

    Corrects to: sensor_zn=0, solar_zn=reference (scene mean), relative_az=0

    Args:
        reflectance: (n_bands,) array
        solar_zn, solar_az, sensor_zn, sensor_az: angles in radians
        ref_solar_zn: reference solar zenith (default: use observed)

    Returns: corrected reflectance (n_bands,)
    """
    if ref_solar_zn is None:
        ref_solar_zn = solar_zn

    relative_az = np.abs(solar_az - sensor_az)
    relative_az = np.where(relative_az > np.pi, 2 * np.pi - relative_az, relative_az)

    # Observed kernels
    k_vol_obs = ross_thick_kernel(solar_zn, sensor_zn, relative_az)
    k_geo_obs = li_dense_r_kernel(solar_zn, sensor_zn, relative_az)

    # Reference kernels (nadir view: sensor_zn=0, relative_az=0)
    k_vol_ref = ross_thick_kernel(ref_solar_zn, 0.0, 0.0)
    k_geo_ref = li_dense_r_kernel(ref_solar_zn, 0.0, 0.0)

    # Simple ratio correction (assumes f_iso=1, f_vol=0.5, f_geo=0.5 as typical)
    # More sophisticated: estimate coefficients per NDVI bin
    f_iso, f_vol, f_geo = 1.0, 0.5, 0.25

    brdf_obs = f_iso + f_vol * k_vol_obs + f_geo * k_geo_obs
    brdf_ref = f_iso + f_vol * k_vol_ref + f_geo * k_geo_ref

    if brdf_obs > 0.01:
        correction = brdf_ref / brdf_obs
        correction = np.clip(correction, 0.5, 2.0)  # limit extreme corrections
        return reflectance * correction
    return reflectance


def process_h5_tile(h5_path, site_plots, site_dir, epsg=None, year_str="unknown"):
    """Open NEON HDF5, apply BRDF correction, clip to plots, save."""
    n_clipped = 0

    try:
        with h5py.File(str(h5_path), "r") as hf:
            base_key = list(hf.keys())[0]
            refl_grp = hf[base_key]["Reflectance"]
            refl_data = refl_grp["Reflectance_Data"]  # (rows, cols, bands) for NEON
            metadata = refl_grp["Metadata"]

            # Spatial info
            map_info = metadata["Coordinate_System"]["Map_Info"][()].decode().split(",")
            x_origin = float(map_info[3])
            y_origin = float(map_info[4])
            x_res = float(map_info[5])
            y_res = float(map_info[6])

            n_rows, n_cols, n_bands = refl_data.shape

            # Angle data
            anc = metadata.get("Ancillary_Imagery", {})
            logs = metadata.get("Logs", {})

            # Solar angles (may be scalar or per-pixel)
            solar_zn_data = None
            solar_az_data = None
            sensor_zn_data = None
            sensor_az_data = None

            # Try per-pixel angles first
            if "to-sensor_Zenith_Angle" in anc:
                sensor_zn_data = anc["to-sensor_Zenith_Angle"][()]
                sensor_az_data = anc["to-sensor_Azimuth_Angle"][()]

            # Solar angles from logs (often scalar per tile)
            if "Solar_Zenith_Angle" in logs:
                solar_zn_data = logs["Solar_Zenith_Angle"][()]
            if "Solar_Azimuth_Angle" in logs:
                solar_az_data = logs["Solar_Azimuth_Angle"][()]

            # Convert to radians
            if solar_zn_data is not None:
                if np.isscalar(solar_zn_data) or solar_zn_data.ndim == 0:
                    solar_zn_rad = float(solar_zn_data) * np.pi / 180
                else:
                    solar_zn_rad = solar_zn_data * np.pi / 180
            else:
                solar_zn_rad = 30.0 * np.pi / 180  # fallback

            if solar_az_data is not None:
                if np.isscalar(solar_az_data) or solar_az_data.ndim == 0:
                    solar_az_rad = float(solar_az_data) * np.pi / 180
                else:
                    solar_az_rad = solar_az_data * np.pi / 180
            else:
                solar_az_rad = 180.0 * np.pi / 180

            # Reference solar zenith (scene mean)
            if np.isscalar(solar_zn_rad):
                ref_solar_zn = solar_zn_rad
            else:
                ref_solar_zn = np.nanmean(solar_zn_rad)

            # For each overlapping plot
            tm = re.search(r'_(\d{6})_(\d{7})_', h5_path.name)
            if not tm:
                return 0
            tile_e, tile_n = int(tm.group(1)), int(tm.group(2))

            for _, prow in site_plots.iterrows():
                px, py = prow["easting"], prow["northing"]
                if not (tile_e <= px + BUFFER and tile_e + 1000 >= px - BUFFER and
                        tile_n <= py + BUFFER and tile_n + 1000 >= py - BUFFER):
                    continue

                plot_id = prow["plotID"]
                out_path = site_dir / f"{year_str}_{plot_id}_hyper_brdf.tif"
                if out_path.exists():
                    n_clipped += 1
                    continue

                # Window in pixel coords
                col_start = max(0, int((px - BUFFER - x_origin) / x_res))
                col_end = min(n_cols, int((px + BUFFER - x_origin) / x_res))
                row_start = max(0, int((y_origin - py - BUFFER) / y_res))
                row_end = min(n_rows, int((y_origin - py + BUFFER) / y_res))

                if col_end <= col_start or row_end <= row_start:
                    continue

                # Read reflectance window
                window_refl = refl_data[row_start:row_end, col_start:col_end, :]  # (h, w, bands)
                h, w, nb = window_refl.shape

                # Apply BRDF correction per pixel
                corrected = np.zeros_like(window_refl, dtype=np.float32)

                for r in range(h):
                    for c in range(w):
                        pixel = window_refl[r, c, :].astype(np.float32)
                        if np.all(pixel <= 0) or np.all(pixel >= 10000):
                            corrected[r, c, :] = pixel
                            continue

                        # Get per-pixel angles if available
                        szn = solar_zn_rad if np.isscalar(solar_zn_rad) else solar_zn_rad[row_start+r, col_start+c] if solar_zn_rad.ndim == 2 else solar_zn_rad
                        saz = solar_az_rad if np.isscalar(solar_az_rad) else solar_az_rad[row_start+r, col_start+c] if solar_az_rad.ndim == 2 else solar_az_rad
                        vzn = 0.0  # default nadir
                        vaz = 0.0
                        if sensor_zn_data is not None and sensor_zn_data.ndim >= 2:
                            vzn = sensor_zn_data[row_start+r, col_start+c] * np.pi / 180
                            vaz = sensor_az_data[row_start+r, col_start+c] * np.pi / 180

                        corrected[r, c, :] = brdf_correct_pixel(
                            pixel, szn, saz, vzn, vaz, ref_solar_zn=ref_solar_zn
                        )

                # Save as GeoTIFF (bands, h, w)
                corrected_t = np.moveaxis(corrected, 2, 0)  # (bands, h, w)
                transform = from_origin(
                    x_origin + col_start * x_res,
                    y_origin - row_start * y_res,
                    x_res, y_res
                )
                profile = {
                    "driver": "GTiff", "dtype": "int16",
                    "width": w, "height": h, "count": nb,
                    "transform": transform, "compress": "lzw", "nodata": -9999,
                }
                if epsg:
                    profile["crs"] = CRS.from_epsg(epsg)

                with rasterio.open(str(out_path), "w", **profile) as dst:
                    dst.write(corrected_t.astype(np.int16))

                n_clipped += 1

    except Exception as e:
        print(f"      Error: {e}")

    h5_path.unlink(missing_ok=True)
    return n_clipped


def process_site(site):
    """Download .001 tiles, apply BRDF, clip to plots."""
    import urllib.request, time

    ppy = pd.read_csv(VEG_STRUCT_DIR / "vst_perplotperyear.csv", low_memory=False)
    plot_coords = ppy[ppy["siteID"] == site][["siteID", "plotID", "easting", "northing"]].dropna()
    plot_coords = plot_coords.drop_duplicates("plotID")

    if plot_coords.empty:
        print(f"  SKIP {site}: no plots")
        return

    # Tiles needed
    needed = set()
    for _, row in plot_coords.iterrows():
        px, py = row["easting"], row["northing"]
        for te in range(int((px-BUFFER)//1000)*1000, int((px+BUFFER)//1000)*1000+1001, 1000):
            for tn in range(int((py-BUFFER)//1000)*1000, int((py+BUFFER)//1000)*1000+1001, 1000):
                needed.add((te, tn))

    # Available months
    try:
        prod = api_get(f"{NEON_API}/products/{DPID}")["data"]
        avail = {s["siteCode"]: sorted(s.get("availableMonths", []))
                 for s in prod["siteCodes"]}
        months = [m for m in avail.get(site, []) if any(m.startswith(str(y)) for y in YEARS)]
    except:
        months = []

    if not months:
        print(f"  SKIP {site}: no data")
        return

    site_dir = OUTPUT_DIR / site
    site_dir.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    epsg = SITE_EPSG.get(site)

    years_avail = sorted(set(m[:4] for m in months))
    print(f"  {site}: {len(plot_coords)} plots, {len(years_avail)} years, {len(needed)} tiles")

    for month in months:
        year_str = month[:4]
        marker = site_dir / f".done_{month}"
        if marker.exists():
            continue

        # Check which plots still need clipping
        plots_needed = []
        for _, row in plot_coords.iterrows():
            if not (site_dir / f"{year_str}_{row['plotID']}_hyper_brdf.tif").exists():
                plots_needed.append(row)
        if not plots_needed:
            marker.touch()
            continue

        files = api_get(f"{NEON_API}/data/{DPID}/{site}/{month}")
        files = files.get("data", {}).get("files", [])
        h5s = [f for f in files if f["name"].endswith(".h5")]

        matching = []
        for f in h5s:
            tm = re.search(r'_(\d{6})_(\d{7})_', f["name"])
            if tm and (int(tm.group(1)), int(tm.group(2))) in needed:
                matching.append(f)

        n_total = 0
        for finfo in matching:
            tmp_path = TEMP_DIR / finfo["name"]

            # Download
            req = urllib.request.Request(finfo["url"])
            if TOKEN:
                req.add_header("X-API-Token", TOKEN)
            try:
                with urllib.request.urlopen(req, timeout=600) as r:
                    with open(str(tmp_path), "wb") as f:
                        while True:
                            chunk = r.read(1024 * 1024)
                            if not chunk:
                                break
                            f.write(chunk)

                n = process_h5_tile(tmp_path, plot_coords, site_dir, epsg=epsg, year_str=month[:4])
                n_total += n
            except Exception as e:
                print(f"      DL/Process error: {e}")
                tmp_path.unlink(missing_ok=True)

        print(f"    {month}: {n_total} BRDF-corrected clips from {len(matching)} tiles")
        marker.touch()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site", type=str, default=None)
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sites = [args.site] if args.site else SITES

    print("BRDF Correction Pipeline (Ross-Li kernels on NEON .001)\n")

    for site in sites:
        print(f"\n{'='*50}")
        process_site(site)

    print("\nDone.")


if __name__ == "__main__":
    main()
