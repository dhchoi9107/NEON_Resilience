"""
Extract MODIS (MOD17A3HGF) + PML-V2 GPP at the 7 NEW-site plots via GEE,
using the IDENTICAL methods as 149_modis_gpp.py / 150_pml_gpp.py, then append
to the existing 19-site GPP files -> _26 versions (existing values untouched).
Usage: python extract_gpp_newsites.py geedankook
"""
import sys, os, pandas as pd, ee
PROJECT = sys.argv[1] if len(sys.argv) > 1 else 'geedankook'
ee.Initialize(project=PROJECT); print("EE init:", PROJECT, flush=True)
D = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
NEW = ["DELA", "LENO", "UKFS", "YELL", "BONA", "DEJU", "HEAL"]

pl = pd.read_csv(os.path.join(D, "plot_lonlat_26.csv"))
pl = pl[pl.siteID.isin(NEW)]
plots = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon, r.lat]), {'plotID': r.plotID})
                              for r in pl.itertuples()])
print("new-site plots:", len(pl), flush=True)

# ── MODIS annual GPP (identical to 149) ──
GPP = ee.ImageCollection("MODIS/061/MOD17A3HGF")
rows = []
for y in range(2016, 2024):
    img = GPP.filterDate(f"{y}-01-01", f"{y}-12-31").first().select('Gpp')
    img = img.updateMask(img.lt(30000)).multiply(0.0001).rename('gpp')
    fc = img.reduceRegions(plots, ee.Reducer.mean(), 500).getInfo()
    for f in fc['features']:
        p = f['properties']; rows.append((p['plotID'], y, p.get('mean', p.get('gpp'))))
    print(f"  MODIS {y} done", flush=True)
gm = pd.DataFrame(rows, columns=['plotID', 'year', 'gpp']).groupby('plotID')['gpp'].mean().reset_index().rename(columns={'gpp': 'modis_gpp'})

# ── PML-V2 GPP (identical to 150) ──
PML = ee.ImageCollection("CAS/IGSNRR/PML/V2_v017").select('GPP')
prows = []
for y in range(2016, 2021):
    img = PML.filterDate(f"{y}-01-01", f"{y}-12-31").mean()
    fc = img.reduceRegions(plots, ee.Reducer.mean(), 500).getInfo()
    for f in fc['features']:
        p = f['properties']; prows.append((p['plotID'], y, p.get('mean', p.get('GPP'))))
    print(f"  PML {y} done", flush=True)
gp = pd.DataFrame(prows, columns=['plotID', 'year', 'pml']).groupby('plotID')['pml'].mean().reset_index().rename(columns={'pml': 'pml_gpp'})

# ── append to existing -> _26 ──
for src, newdf, col in [("plot_modis_gpp.csv", gm, "modis_gpp"), ("plot_pml_gpp.csv", gp, "pml_gpp")]:
    exist = pd.read_csv(os.path.join(D, src))
    comb = pd.concat([exist, newdf], ignore_index=True).drop_duplicates('plotID')
    out = os.path.join(D, src.replace(".csv", "_26.csv"))
    comb.to_csv(out, index=False)
    print(f"-> {out}: {len(comb)} plots, new-site non-null {col}: {newdf[col].notna().sum()}/{len(newdf)}", flush=True)

print("\nnew-site GPP medians by site:", flush=True)
mm = gm.merge(gp, on='plotID', how='outer'); mm['siteID'] = mm.plotID.str.split('_').str[0]
print(mm.groupby('siteID').agg(modis=('modis_gpp', 'median'), pml=('pml_gpp', 'median')).round(1).to_string(), flush=True)
