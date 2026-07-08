"""
NEON_v2 STEP 1 — Re-extract plot VIs DIRECTLY from NEON .002 bidirectional
(BRDF + topo corrected) VI products. Authoritative, single-source.

Only years with the NEON .002 BRDF VI product exist on disk
(2016,2017,2018,2022,2023,2024,2025,2026) -> this script can ONLY produce
those years, so the output is guaranteed 100% NEON-BRDF source. No hyperspectral
fallback, no mixing. Every row tagged source='NEON_DP3.30026.002_bidirectional'.

Output: NEON_v2/data/plot_vi_neon_brdf.csv  (plot-level mean VI at 10 m)
Does NOT overwrite any original file.
"""
import sys, os, re, io, zipfile
sys.path.insert(0, r"C:\Users\star1\Documents\GitHub\NEON_Resilience")
import numpy as np, pandas as pd
from compute.compute_plot_spectral_1m import (
    discover_tiles_by_year, find_tile_for_point, extract_plot_1m,
    aggregate_to_grain, VI_BANDS,
)
from site_config import SITES, VEG_STRUCT_DIR

KEEP = ["NDVI", "EVI", "ARVI", "SAVI"]   # PRI/fPAR excluded (cross-year artifact)
GRAIN = 10
SOURCE = "NEON_DP3.30026.002_bidirectional"
OUT = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_vi_neon_brdf.csv"

# plot coordinates (one per plot)
ppy = pd.read_csv(VEG_STRUCT_DIR / "vst_perplotperyear.csv", low_memory=False)
coords = ppy[["siteID", "plotID", "easting", "northing"]].dropna().drop_duplicates("plotID")
print(f"plots with coords: {len(coords)}", flush=True)

rows = []
for site in SITES:
    sc = coords[coords["siteID"] == site]
    if sc.empty:
        continue
    tiles = discover_tiles_by_year(site)   # only finds .002 zips on disk
    if not tiles:
        continue
    yrs = sorted(tiles.keys())
    print(f"{site}: years {yrs}  ({len(sc)} plots)", flush=True)
    for year in yrs:
        yt = tiles[year]
        for _, pr in sc.iterrows():
            x, y = float(pr.easting), float(pr.northing)
            entry = find_tile_for_point(site, x, y, yt)
            if entry is None or "VI_zip" not in entry:
                continue
            stack = extract_plot_1m(entry, x, y)   # (n_bands, H, W), bands = VI_BANDS+[LAI,fPAR]
            if stack is None:
                continue
            cells = aggregate_to_grain(stack, GRAIN)   # (n_cells, n_bands)
            if len(cells) == 0:
                continue
            rec = {"siteID": site, "plotID": pr.plotID, "year": year,
                   "grain_m": GRAIN, "n_cells": len(cells), "source": SOURCE}
            with np.errstate(all="ignore"):
                mean = np.nanmean(cells, axis=0); sd = np.nanstd(cells, axis=0)
            for b in KEEP:
                bi = VI_BANDS.index(b)
                rec[f"{b}_mean"] = float(mean[bi]); rec[f"{b}_sd"] = float(sd[bi])
            rows.append(rec)
        print(f"  {site} {year}: {sum(1 for r in rows if r['siteID']==site and r['year']==year)} plots", flush=True)

out = pd.DataFrame(rows)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
out.to_csv(OUT, index=False)
print(f"\nDONE -> {OUT}", flush=True)
print(f"rows: {len(out)} | years: {sorted(out['year'].unique())} | plots: {out['plotID'].nunique()}", flush=True)
print(out.groupby('year').size().to_string(), flush=True)
