"""
Rebuild plot_level_complete.csv from ALL 1m FSD rasters
========================================================
Reads each plot's 140x140 FSD_1m.tif, extracts 40x40m center,
computes mean/sd per band, merges with pooled taxonomic + LCBD.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import pandas as pd
import rasterio
import warnings
warnings.filterwarnings("ignore")

from site_config import FSD_BANDS

FSD_DIR = Path("E:/neon_lidar/structural_diversity_1m_plots")
OUT_PATH = Path("E:/neon_lidar/model_results/plot_level_complete.csv")

# Band names for columns
band_names = [
    "rumple", "top_rugosity", "mean_max_canopy_ht", "max_canopy_ht",
    "deepgap_fraction", "meanH", "vert_sd", "vertCV",
    "mean_sd", "sd_sd", "GC", "GFP", "VCI",
    "q25", "q50", "q75", "q95", "HeightRatio",
    "FHD", "LAI", "LAI_subcanopy"
]

WINDOW_SIZE = 140
PLOT_HALF = 20  # 40m plot = ±20m from center
CENTER_START = (WINDOW_SIZE // 2) - PLOT_HALF  # 50
CENTER_END = (WINDOW_SIZE // 2) + PLOT_HALF    # 90

rows = []
n_files = 0
n_errors = 0

for site_dir in sorted(FSD_DIR.iterdir()):
    if not site_dir.is_dir():
        continue
    site = site_dir.name
    for tif in sorted(site_dir.glob("*_FSD_1m.tif")):
        # Parse year and plotID from filename: {year}_{plotID}_FSD_1m.tif
        parts = tif.stem.replace("_FSD_1m", "").split("_", 1)
        if len(parts) != 2:
            continue
        year = int(parts[0])
        plot_id = parts[1]

        try:
            with rasterio.open(str(tif)) as src:
                data = src.read()  # (21, 140, 140)
        except Exception:
            n_errors += 1
            continue

        # Extract 40×40m center
        center = data[:, CENTER_START:CENTER_END, CENTER_START:CENTER_END]  # (21, 40, 40)

        row = {"siteID": site, "plotID": plot_id, "fsd_year": year}
        for i, bname in enumerate(band_names):
            band = center[i]
            valid = band[~np.isnan(band)]
            if len(valid) > 0:
                row[f"{bname}_mean"] = float(np.nanmean(valid))
                row[f"{bname}_sd"] = float(np.nanstd(valid))
            else:
                row[f"{bname}_mean"] = np.nan
                row[f"{bname}_sd"] = np.nan

        rows.append(row)
        n_files += 1

    print(f"  {site}: {sum(1 for r in rows if r['siteID'] == site)} plot-years")

fsd_df = pd.DataFrame(rows)
print(f"\nTotal: {n_files} plot-years, {n_errors} errors")
print(f"Years: {sorted(fsd_df.fsd_year.unique())}")
print(f"Sites: {fsd_df.siteID.nunique()}")

# Merge with pooled taxonomic
alpha = pd.read_csv("E:/neon_lidar/taxonomic_diversity/alpha_diversity_pooled.csv")
alpha = alpha[["siteID", "plotID", "richness", "shannon", "simpson"]]

# Merge with LCBD
lcbd_src = pd.read_csv("E:/neon_lidar/model_results/plot_level_dbh10.csv",
                        usecols=["siteID", "plotID", "lcbd_bray", "lcbd_turnover", "lcbd_nestedness"])
lcbd_src = lcbd_src.drop_duplicates("plotID")

merged = fsd_df.merge(alpha, on=["siteID", "plotID"], how="left")
merged = merged.merge(lcbd_src, on=["siteID", "plotID"], how="left")

# Backup old file
if OUT_PATH.exists():
    backup = OUT_PATH.with_suffix(".csv.bak")
    OUT_PATH.rename(backup)
    print(f"\nBacked up old file to {backup.name}")

merged.to_csv(str(OUT_PATH), index=False)
print(f"Saved: {OUT_PATH} ({len(merged)} rows)")
print(f"Year distribution:")
print(merged.groupby("fsd_year").size())
print("\nDone.")
