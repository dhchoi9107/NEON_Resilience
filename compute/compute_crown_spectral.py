"""
Crown-level Spectral Diversity
================================
Delineates individual tree crowns from 1m CHM (watershed),
extracts mean hyperspectral signature per crown,
computes crown-level spectral diversity metrics.

Hypothesis: crown-level diversity removes within-crown noise
and should better predict species diversity than pixel-level.

Usage:
  python compute/compute_crown_spectral.py
  python compute/compute_crown_spectral.py --site HARV
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import warnings
import re

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
from scipy.ndimage import label, maximum_filter
from skimage.segmentation import watershed
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

from site_config import SITES, SPEC_DIV_DIR
from compute.compute_hyperspectral_diversity import get_good_band_mask

warnings.filterwarnings("ignore")

FSD_1M_DIR = Path("E:/neon_lidar/structural_diversity_1m_plots")
HYPER_DIR = Path("E:/neon_lidar/hyperspectral_plots")
MIN_CROWN_PX = 3
MIN_TREE_HT = 2.0


def delineate_crowns(chm, window_size=5):
    """Watershed crown segmentation from CHM."""
    chm = chm.copy()
    chm[~np.isfinite(chm)] = 0
    chm[chm < MIN_TREE_HT] = 0

    local_max = maximum_filter(chm, footprint=np.ones((window_size, window_size)))
    markers, n_tops = label((chm == local_max) & (chm > MIN_TREE_HT))
    crowns = watershed(-chm, markers, mask=(chm > MIN_TREE_HT))
    return crowns, n_tops


def extract_crown_spectra(crowns, hyper_good, valid_pixel):
    """Extract mean spectrum per crown. Returns (n_crowns, n_bands) array."""
    spectra = []
    for cid in np.unique(crowns):
        if cid == 0:
            continue
        mask = (crowns == cid) & valid_pixel
        if mask.sum() < MIN_CROWN_PX:
            continue
        spectra.append(hyper_good[:, mask].mean(axis=1))
    return np.array(spectra) if spectra else None


def compute_diversity(data, label=""):
    """Compute diversity metrics from (n_samples, n_features) array."""
    if data is None or len(data) < 5:
        return {f"{label}n": 0, f"{label}FRic": np.nan,
                f"{label}RaoQ": np.nan, f"{label}FDiv": np.nan}

    std = StandardScaler().fit_transform(data)
    pca = PCA(n_components=min(5, std.shape[1]))
    scores = pca.fit_transform(std)

    result = {f"{label}n": len(data)}

    # FRic
    try:
        hull = ConvexHull(scores[:min(len(scores), 2000), :3])
        result[f"{label}FRic"] = float(hull.volume)
    except:
        result[f"{label}FRic"] = np.nan

    # Rao's Q
    sub = scores[:min(len(scores), 500)]
    result[f"{label}RaoQ"] = float(np.mean(pdist(sub)))

    # FDiv
    centroid = np.mean(scores, axis=0)
    dists = np.sqrt(np.sum((scores - centroid) ** 2, axis=1))
    result[f"{label}FDiv"] = float(np.mean(dists) / max(np.max(dists), 1e-8))

    return result


def process_plot(fsd_path, hyper_path):
    """Process one plot: crown delineation + crown/pixel diversity."""
    # CHM for crown delineation
    with rasterio.open(str(fsd_path)) as src:
        chm = src.read(3).astype(np.float32)  # mean_max_canopy_ht
        fsd_tf, fsd_b = src.transform, src.bounds

    crowns, n_tops = delineate_crowns(chm)

    # Hyperspectral
    with rasterio.open(str(hyper_path)) as src:
        hyper = src.read().astype(np.float32)
        hyper_tf, hyper_b = src.transform, src.bounds
        n_bands = src.count

    # Overlap region
    ox0 = max(fsd_b.left, hyper_b.left)
    ox1 = min(fsd_b.right, hyper_b.right)
    oy0 = max(fsd_b.bottom, hyper_b.bottom)
    oy1 = min(fsd_b.top, hyper_b.top)

    if ox1 <= ox0 or oy1 <= oy0:
        return None

    fr1, fc1 = rowcol(fsd_tf, ox0, oy1)
    fr2, fc2 = rowcol(fsd_tf, ox1, oy0)
    hr1, hc1 = rowcol(hyper_tf, ox0, oy1)
    hr2, hc2 = rowcol(hyper_tf, ox1, oy0)
    h = min(fr2 - fr1, hr2 - hr1)
    w = min(fc2 - fc1, hc2 - hc1)

    if h < 10 or w < 10:
        return None

    crowns_c = crowns[fr1:fr1+h, fc1:fc1+w]
    hyper_c = hyper[:, hr1:hr1+h, hc1:hc1+w]

    good_mask, _ = get_good_band_mask(n_bands)
    hyper_g = hyper_c[good_mask]
    valid_px = np.all((hyper_g > 0) & (hyper_g < 10000), axis=0)

    # Crown-level
    crown_spectra = extract_crown_spectra(crowns_c, hyper_g, valid_px)
    crown_div = compute_diversity(crown_spectra, label="crown_")

    # Pixel-level (same area for fair comparison)
    px_flat = hyper_g.reshape(hyper_g.shape[0], -1).T
    px_valid = px_flat[valid_px.ravel()]
    if len(px_valid) > 2000:
        idx = np.random.choice(len(px_valid), 2000, replace=False)
        px_valid = px_valid[idx]
    pixel_div = compute_diversity(px_valid, label="pixel_")

    return {**crown_div, **pixel_div, "n_tree_tops": n_tops}


def process_site(site):
    """Process all plots for a site."""
    hyper_dir = HYPER_DIR / site
    fsd_dir = FSD_1M_DIR / site

    if not hyper_dir.exists() or not fsd_dir.exists():
        print(f"  SKIP {site}")
        return []

    hyper_files = {f.stem.replace("_hyper", ""): f
                   for f in hyper_dir.glob("*_hyper.tif")}
    results = []

    for key, hyper_path in sorted(hyper_files.items()):
        m = re.match(r"(\d{4})_(.+)", key)
        if not m:
            continue
        year, plot_id = m.group(1), m.group(2)

        fsd_path = fsd_dir / f"{year}_{plot_id}_FSD_1m.tif"
        if not fsd_path.exists():
            # Try any year
            candidates = sorted(fsd_dir.glob(f"*_{plot_id}_FSD_1m.tif"), reverse=True)
            if not candidates:
                continue
            fsd_path = candidates[0]

        row = process_plot(fsd_path, hyper_path)
        if row is None:
            continue

        row["siteID"] = site
        row["plotID"] = plot_id
        row["year"] = int(year)
        results.append(row)

    if results:
        mean_cn = np.mean([r["crown_n"] for r in results])
        mean_cf = np.nanmean([r["crown_FRic"] for r in results])
        print(f"  {site}: {len(results)} plots, "
              f"mean {mean_cn:.0f} crowns/plot, crown FRic={mean_cf:.0f}")

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site", type=str, default=None)
    args = parser.parse_args()

    sites = [args.site] if args.site else SITES
    all_results = []

    print("Crown-level vs Pixel-level Spectral Diversity\n")

    for site in sites:
        rows = process_site(site)
        all_results.extend(rows)

    if all_results:
        df = pd.DataFrame(all_results)
        out_path = SPEC_DIV_DIR / "crown_spectral_diversity.csv"
        df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path} ({len(df)} plots)")

        # Quick comparison
        print("\n=== Crown vs Pixel (mean across all plots) ===")
        for metric in ["FRic", "RaoQ", "FDiv"]:
            cm = df[f"crown_{metric}"].mean()
            pm = df[f"pixel_{metric}"].mean()
            print(f"  {metric}: crown={cm:.1f}, pixel={pm:.1f}, ratio={cm/pm:.2f}")


if __name__ == "__main__":
    main()
