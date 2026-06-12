"""
DeepForest Crown Detection + Crown-level Spectral Diversity
=============================================================
Downloads NEON RGB tiles (DP3.30010.001) for plot locations,
runs DeepForest pre-trained model for crown detection,
then extracts mean hyperspectral signature per detected crown.

Pipeline:
  1. Download RGB tile → clip to plot ± 50m → delete tile
  2. Run DeepForest → bounding boxes per crown
  3. Match crowns with hyperspectral data → mean spectrum per crown
  4. Compute crown-level spectral diversity

Usage:
  python compute/compute_deepforest_crowns.py --site HARV
  python compute/compute_deepforest_crowns.py  # all sites
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import re
import warnings

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol, from_origin
from rasterio.crs import CRS
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

from site_config import (SITES, YEARS, VEG_STRUCT_DIR, TOKEN, NEON_API,
                          api_get, SITE_EPSG)
from compute.compute_hyperspectral_diversity import get_good_band_mask

warnings.filterwarnings("ignore")

DPID_RGB = "DP3.30010.001"
BUFFER = 50
RGB_DIR = Path("E:/neon_lidar/rgb_plots")
HYPER_DIR = Path("E:/neon_lidar/hyperspectral_plots")
OUTPUT_DIR = Path("E:/neon_lidar/spectral_diversity")


def download_and_clip_rgb(site, plot_coords, year_month):
    """Download RGB tiles, clip to plot windows, delete originals."""
    import urllib.request, time

    site_dir = RGB_DIR / site
    site_dir.mkdir(parents=True, exist_ok=True)
    year_str = year_month[:4]

    # Check what clips already exist
    needed_tiles = set()
    plot_tile_map = {}  # tile_key -> [plot rows]
    for _, row in plot_coords.iterrows():
        px, py = row["easting"], row["northing"]
        out_path = site_dir / f"{year_str}_{row['plotID']}_rgb.tif"
        if out_path.exists():
            continue
        tile_e = int(px // 1000) * 1000
        tile_n = int(py // 1000) * 1000
        key = (tile_e, tile_n)
        needed_tiles.add(key)
        if key not in plot_tile_map:
            plot_tile_map[key] = []
        plot_tile_map[key].append(row)

    if not needed_tiles:
        return

    # Get file list
    files = api_get(f"{NEON_API}/data/{DPID_RGB}/{site}/{year_month}")
    files = files.get("data", {}).get("files", [])
    tifs = [f for f in files if f["name"].endswith(".tif") and "image" in f["name"].lower()]

    for finfo in tifs:
        m = re.search(r'_(\d{6})_(\d{7})_', finfo["name"])
        if not m:
            continue
        te, tn = int(m.group(1)), int(m.group(2))
        key = (te, tn)
        if key not in plot_tile_map:
            continue

        # Download
        tmp_path = site_dir / f"_tmp_{finfo['name']}"
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
        except Exception as e:
            print(f"      DL fail: {e}")
            tmp_path.unlink(missing_ok=True)
            continue

        # Clip to each plot
        try:
            with rasterio.open(str(tmp_path)) as src:
                for prow in plot_tile_map[key]:
                    px, py = prow["easting"], prow["northing"]
                    out_path = site_dir / f"{year_str}_{prow['plotID']}_rgb.tif"
                    if out_path.exists():
                        continue

                    window = rasterio.windows.from_bounds(
                        px - BUFFER, py - BUFFER, px + BUFFER, py + BUFFER,
                        transform=src.transform
                    )
                    window = window.intersection(
                        rasterio.windows.Window(0, 0, src.width, src.height)
                    )
                    if window.width < 10 or window.height < 10:
                        continue

                    data = src.read(window=window)
                    transform = rasterio.windows.transform(window, src.transform)
                    profile = src.profile.copy()
                    profile.update(width=int(window.width), height=int(window.height),
                                   transform=transform, compress="lzw")

                    with rasterio.open(str(out_path), "w", **profile) as dst:
                        dst.write(data)
        except Exception as e:
            print(f"      Clip fail: {e}")

        tmp_path.unlink(missing_ok=True)


_DF_MODEL = None

def get_deepforest_model():
    global _DF_MODEL
    if _DF_MODEL is None:
        from deepforest import main as df_main
        _DF_MODEL = df_main.deepforest()
        _DF_MODEL.load_model("weecology/deepforest-tree")
    return _DF_MODEL


def run_deepforest_on_plot(rgb_path):
    """Run DeepForest tree detection on a plot RGB clip.

    Returns DataFrame with xmin, ymin, xmax, ymax, score columns (pixel coords).
    """
    m = get_deepforest_model()

    boxes = m.predict_tile(
        path=str(rgb_path),
        patch_size=400,
        patch_overlap=0.15,
    )

    if boxes is None or len(boxes) == 0:
        return pd.DataFrame()

    return boxes


def crowns_to_spectral(boxes, rgb_path, hyper_path):
    """Extract mean hyperspectral signature per DeepForest crown.

    Converts RGB pixel bounding boxes to UTM coordinates,
    matches with hyperspectral data.
    Returns (n_crowns, n_bands) array.
    """
    with rasterio.open(str(rgb_path)) as src:
        rgb_tf = src.transform
        rgb_res = src.res[0]  # typically 0.1m

    with rasterio.open(str(hyper_path)) as src:
        hyper = src.read().astype(np.float32)
        hyper_tf = src.transform
        hyper_res = src.res[0]  # 1m
        n_bands = src.count

    good_mask, _ = get_good_band_mask(n_bands)
    hyper_g = hyper[good_mask]
    valid_px = np.all((hyper_g > 0) & (hyper_g < 10000), axis=0)

    crown_spectra = []
    for _, box in boxes.iterrows():
        # Convert RGB pixel coords to UTM
        # RGB is 0.1m, so pixel (col, row) → UTM
        cx = rgb_tf.c + (box["xmin"] + box["xmax"]) / 2 * rgb_res
        cy = rgb_tf.f - (box["ymin"] + box["ymax"]) / 2 * rgb_res

        # Crown radius estimate (half of bbox diagonal in meters)
        bw = (box["xmax"] - box["xmin"]) * rgb_res
        bh = (box["ymax"] - box["ymin"]) * rgb_res
        radius = max(bw, bh) / 2

        # Find corresponding pixels in hyperspectral (1m)
        hcol = int((cx - hyper_tf.c) / hyper_res)
        hrow = int((hyper_tf.f - cy) / hyper_res)
        hr = max(1, int(radius / hyper_res))

        r1 = max(0, hrow - hr)
        r2 = min(hyper_g.shape[1], hrow + hr + 1)
        c1 = max(0, hcol - hr)
        c2 = min(hyper_g.shape[2], hcol + hr + 1)

        if r2 <= r1 or c2 <= c1:
            continue

        patch = hyper_g[:, r1:r2, c1:c2]
        patch_valid = valid_px[r1:r2, c1:c2]

        if patch_valid.sum() < 1:
            continue

        mean_spec = patch[:, patch_valid].mean(axis=1)
        crown_spectra.append(mean_spec)

    return np.array(crown_spectra) if crown_spectra else None


def compute_diversity(data, label=""):
    if data is None or len(data) < 5:
        return {f"{label}n": len(data) if data is not None else 0,
                f"{label}FRic": np.nan, f"{label}RaoQ": np.nan, f"{label}FDiv": np.nan}

    std = StandardScaler().fit_transform(data)
    pca = PCA(n_components=min(5, std.shape[1]))
    scores = pca.fit_transform(std)

    result = {f"{label}n": len(data)}
    try:
        hull = ConvexHull(scores[:min(2000, len(scores)), :3])
        result[f"{label}FRic"] = float(hull.volume)
    except:
        result[f"{label}FRic"] = np.nan

    sub = scores[:min(500, len(scores))]
    result[f"{label}RaoQ"] = float(np.mean(pdist(sub)))

    centroid = np.mean(scores, axis=0)
    dists = np.sqrt(np.sum((scores - centroid) ** 2, axis=1))
    result[f"{label}FDiv"] = float(np.mean(dists) / max(np.max(dists), 1e-8))

    return result


def process_site(site):
    """Full pipeline: download RGB → DeepForest → crown spectra → diversity."""
    ppy = pd.read_csv(VEG_STRUCT_DIR / "vst_perplotperyear.csv", low_memory=False)
    plot_coords = ppy[ppy["siteID"] == site][["siteID", "plotID", "easting", "northing"]].dropna()
    plot_coords = plot_coords.drop_duplicates("plotID")

    if plot_coords.empty:
        print(f"  SKIP {site}: no plots")
        return []

    hyper_dir = HYPER_DIR / site
    if not hyper_dir.exists():
        print(f"  SKIP {site}: no hyperspectral")
        return []

    # Find available RGB months (latest year)
    try:
        prod = api_get(f"{NEON_API}/products/{DPID_RGB}")["data"]
        avail = {}
        for s in prod["siteCodes"]:
            avail[s["siteCode"]] = sorted(s.get("availableMonths", []))
        months = [m for m in avail.get(site, []) if any(m.startswith(str(y)) for y in YEARS)]
    except:
        months = []

    if not months:
        print(f"  SKIP {site}: no RGB data")
        return []

    # All available years (not just latest)
    rgb_months = months
    years_avail = sorted(set(m[:4] for m in rgb_months))

    print(f"  {site}: {len(plot_coords)} plots, RGB {years_avail}")

    # Step 1: Download + clip RGB (all years)
    print(f"    Downloading RGB...")
    for month in rgb_months:
        download_and_clip_rgb(site, plot_coords, month)

    rgb_dir = RGB_DIR / site
    rgb_files = sorted(rgb_dir.glob("*_rgb.tif"))

    # Step 2-3: DeepForest + crown spectra (per year per plot)
    results = []
    for rgb_path in rgb_files:
        m = re.match(r"(\d{4})_(.+)_rgb\.tif", rgb_path.name)
        if not m:
            continue
        year, plot_id = m.group(1), m.group(2)

        # DeepForest detection (RGB only — no hyperspectral needed)
        boxes = run_deepforest_on_plot(rgb_path)
        if boxes.empty:
            continue

        row = {
            "siteID": site,
            "plotID": plot_id,
            "year": int(year),
            "n_detections": len(boxes),
            "mean_score": float(boxes["score"].mean()) if "score" in boxes.columns else np.nan,
            "mean_crown_area_m2": float(
                ((boxes["xmax"] - boxes["xmin"]) * (boxes["ymax"] - boxes["ymin"])).mean()
                * abs(rasterio.open(str(rgb_path)).transform.a) ** 2
            ) if len(boxes) > 0 else np.nan,
        }

        # Crown spectra (only if hyperspectral available)
        hyper_path = hyper_dir / f"{year}_{plot_id}_hyper.tif"
        if not hyper_path.exists():
            candidates = sorted(hyper_dir.glob(f"*_{plot_id}_hyper.tif"))
            hyper_path = candidates[0] if candidates else None

        if hyper_path and hyper_path.exists():
            crown_spectra = crowns_to_spectral(boxes, rgb_path, hyper_path)
            div = compute_diversity(crown_spectra, label="df_crown_")
            row.update(div)
        else:
            row["df_crown_n"] = 0

        results.append(row)

    if results:
        n_years = len(set(r["year"] for r in results))
        mean_n = np.mean([r["n_detections"] for r in results])
        print(f"    {len(results)} plot-years ({n_years} years), {mean_n:.0f} detections/plot")

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site", type=str, default=None)
    args = parser.parse_args()

    sites = [args.site] if args.site else SITES
    all_results = []

    print("DeepForest Crown Detection + Spectral Diversity\n")

    for site in sites:
        print(f"\n{'='*50}")
        rows = process_site(site)
        all_results.extend(rows)

    if all_results:
        df = pd.DataFrame(all_results)
        out_path = OUTPUT_DIR / "deepforest_crown_diversity.csv"
        df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path} ({len(df)} plots)")


if __name__ == "__main__":
    main()
