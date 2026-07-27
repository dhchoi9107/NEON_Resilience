"""
Extract NEON .002 BRDF VI for TREE ONLY and MERGE into plot_vi_neon_brdf.csv.
Mirrors NEON_v2/scripts/01_extract_vi_neon_brdf.py extraction logic, but:
  - scoped to TREE (uses SITES_ENV=TREE via site_config override),
  - does NOT overwrite the shared file: backs up, drops any existing TREE rows,
    appends freshly extracted TREE rows, writes merged result.
Safe to re-run (idempotent for TREE).
"""
import sys, os, shutil
sys.path.insert(0, r"C:\Users\star1\Documents\GitHub\NEON_Resilience")
import numpy as np, pandas as pd
from compute.compute_plot_spectral_1m import (
    discover_tiles_by_year, find_tile_for_point, extract_plot_1m,
    aggregate_to_grain, VI_BANDS,
)
from site_config import VEG_STRUCT_DIR

SITE = "TREE"
# TREE (Treehaven) shares its NEON AOP flight box with STEI (Steigerwaldt); the
# airborne tiles are published under the STEI site code (NEON_D05_STEI_*_VegIndices.zip).
# TREE plot coordinates fall inside these STEI tiles, so we discover tiles for STEI
# but keep siteID=TREE. find_tile_for_point() matches purely on UTM coords (site arg unused).
TILE_SITE = "STEI"
KEEP = ["NDVI", "EVI", "ARVI", "SAVI"]
GRAIN = 10
SOURCE = "NEON_DP3.30026.002_bidirectional"
OUT = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_vi_neon_brdf.csv"

ppy = pd.read_csv(VEG_STRUCT_DIR / "vst_perplotperyear.csv", low_memory=False)
coords = ppy[["siteID", "plotID", "easting", "northing"]].dropna().drop_duplicates("plotID")
sc = coords[coords["siteID"] == SITE]
print(f"{SITE}: {len(sc)} plots with coords", flush=True)

tiles = discover_tiles_by_year(TILE_SITE)   # STEI flight box covers TREE plots
if not tiles:
    print(f"ERROR: no .002 VI tiles on disk for flight-box {TILE_SITE} — aborting.", flush=True)
    sys.exit(2)
yrs = sorted(tiles.keys())
print(f"{SITE} via {TILE_SITE} flight box: years on disk = {yrs}", flush=True)

rows = []
for year in yrs:
    yt = tiles[year]
    for _, pr in sc.iterrows():
        x, y = float(pr.easting), float(pr.northing)
        entry = find_tile_for_point(SITE, x, y, yt)
        if entry is None or "VI_zip" not in entry:
            continue
        try:
            stack = extract_plot_1m(entry, x, y)
        except Exception as e:
            print(f"  [skip] {pr.plotID} {year}: {type(e).__name__} ({e})", flush=True)
            continue
        if stack is None:
            continue
        cells = aggregate_to_grain(stack[:len(VI_BANDS)], GRAIN)   # VI bands only (no LAI/fPAR tiles needed)
        if len(cells) == 0:
            continue
        rec = {"siteID": SITE, "plotID": pr.plotID, "year": year,
               "grain_m": GRAIN, "n_cells": len(cells), "source": SOURCE}
        with np.errstate(all="ignore"):
            mean = np.nanmean(cells, axis=0); sd = np.nanstd(cells, axis=0)
        for b in KEEP:
            bi = VI_BANDS.index(b)
            rec[f"{b}_mean"] = float(mean[bi]); rec[f"{b}_sd"] = float(sd[bi])
        rows.append(rec)
    print(f"  {SITE} {year}: {sum(1 for r in rows if r['year']==year)} plots", flush=True)

new = pd.DataFrame(rows)
if new.empty:
    print(f"WARNING: 0 {SITE} VI rows extracted — leaving {OUT} unchanged.", flush=True)
    sys.exit(3)

# ---- MERGE (never blind-overwrite) ----
existing = pd.read_csv(OUT)
shutil.copy2(OUT, OUT + ".bak")
print(f"backup -> {OUT}.bak ({len(existing)} rows, {existing.siteID.nunique()} sites)", flush=True)

kept = existing[existing["siteID"] != SITE]            # drop any stale TREE rows (expect 0)
merged = pd.concat([kept, new], ignore_index=True)
merged = merged.sort_values(["siteID", "plotID", "year"]).reset_index(drop=True)
merged = merged[existing.columns]                       # preserve column order
merged.to_csv(OUT, index=False)

print(f"\nDONE -> {OUT}", flush=True)
print(f"  before: {len(existing)} rows, {existing.siteID.nunique()} sites", flush=True)
print(f"  {SITE} added: {len(new)} rows, {new.plotID.nunique()} plots, years {sorted(new.year.unique())}", flush=True)
print(f"  after:  {len(merged)} rows, {merged.siteID.nunique()} sites", flush=True)
