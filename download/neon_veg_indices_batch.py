#!/usr/bin/env python3
"""
NEON AOP Vegetation Indices & Related Products Batch Downloader
================================================================
Downloads all vegetation-related spectral products (Level 3 mosaic, 1km tiles)
for the same sites/years as the LiDAR structural diversity analysis.

Products:
  DP3.30026.001 - Vegetation Indices (NDVI, EVI, ARVI, PRI, NDLI, NDNI, SAVI)
  DP3.30012.001 - LAI (Leaf Area Index)
  DP3.30014.001 - fPAR (Fraction of Photosynthetically Active Radiation)
  DP3.30011.001 - Albedo
  DP3.30019.001 - Canopy Water Indices
  DP3.30015.001 - CHM (Canopy Height Model)

Usage:
  python neon_veg_indices_batch.py estimate              # check availability & size
  python neon_veg_indices_batch.py download               # download all
  python neon_veg_indices_batch.py download --product VI  # download only Veg Indices
"""
import os, sys, json, time, argparse, urllib.request, urllib.error

# ─── Products ────────────────────────────────────────────────
PRODUCTS = {
    # .002 = BRDF + topographic corrected (bidirectional reflectance)
    "VI":    ("DP3.30026.002", "Vegetation Indices - BRDF corrected"),
    "LAI":   ("DP3.30012.002", "LAI - BRDF corrected"),
    "fPAR":  ("DP3.30014.002", "fPAR - BRDF corrected"),
    # ALB .002 excluded — only 5 site-months available
    "CWI":   ("DP3.30019.002", "Canopy Water Indices - BRDF corrected"),
    # CHM excluded per user request
}

# ─── Sites & Years (same as LiDAR analysis) ─────────────────
SITES = sorted([
    "HARV","BART","UNDE","BLAN","SCBI","SERC","MLBS","ORNL","GRSM","WREF","ABBY",
    "SOAP","TEAK","RMNP","JERC","OSBS","TALL","STEI","TREE",
    "CHEQ",  # included in LiDAR processing
])
YEARS = list(range(2013, 2027))

SAVEPATH = os.environ.get("SAVEPATH", "E:/neon_lidar/vegetation_indices")

# Token
TOKEN = os.environ.get("NEON_TOKEN")
if not TOKEN:
    for _tf in (os.path.join(os.path.dirname(os.path.abspath(__file__)), ".neon_token"),
                os.path.expanduser("~/.neon_token")):
        if os.path.exists(_tf):
            TOKEN = open(_tf).read().strip()
            break

if os.environ.get("SITES_ENV"):
    SITES = [s.strip() for s in os.environ["SITES_ENV"].split(",") if s.strip()]

API = "https://data.neonscience.org/api/v0"

def _get(url):
    req = urllib.request.Request(url)
    if TOKEN:
        req.add_header("X-API-Token", TOKEN)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < 2:
                time.sleep(5)
            else:
                raise

_PROD_CACHE = {}
def available_months(dpid, site):
    if dpid not in _PROD_CACHE:
        d = _get(f"{API}/products/{dpid}")["data"]
        _PROD_CACHE[dpid] = {}
        for s in d["siteCodes"]:
            _PROD_CACHE[dpid][s["siteCode"]] = sorted(s.get("availableMonths", []))
    return _PROD_CACHE[dpid].get(site, [])

def target_months(dpid, site):
    return [ym for ym in available_months(dpid, site) if int(ym[:4]) in YEARS]


# ── ESTIMATE ─────────────────────────────────────────────────
def estimate(product_filter=None):
    products = PRODUCTS if not product_filter else {k: v for k, v in PRODUCTS.items() if k == product_filter}

    for key, (dpid, desc) in products.items():
        print(f"\n{'='*60}")
        print(f"[{key}] {dpid} - {desc}")
        print(f"{'='*60}")
        print(f"{'site':6} {'available months':30} {'count':>6}")

        total_months = 0
        available_sites = 0
        for site in SITES:
            yms = target_months(dpid, site)
            if yms:
                available_sites += 1
                total_months += len(yms)
                print(f"{site:6} {','.join(yms):30} {len(yms):>6}")
            else:
                print(f"{site:6} {'(not available)':30}")
            time.sleep(0.05)

        print(f"  -> {available_sites} sites, {total_months} site-months available")

    print(f"\nSave path: {SAVEPATH}")


# ── DOWNLOAD ─────────────────────────────────────────────────
def download(product_filter=None):
    import neonutilities as nu

    products = PRODUCTS if not product_filter else {k: v for k, v in PRODUCTS.items() if k == product_filter}
    os.makedirs(SAVEPATH, exist_ok=True)
    log = open(os.path.join(SAVEPATH, "download_log.txt"), "a")

    def L(m):
        print(m); log.write(m + "\n"); log.flush()

    L(f"\n{'='*60}")
    L(f"Download started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    L(f"Products: {', '.join(products.keys())}")
    L(f"Sites: {', '.join(SITES)}")
    L(f"Save path: {SAVEPATH}")
    L(f"{'='*60}\n")

    for key, (dpid, desc) in products.items():
        L(f"\n--- [{key}] {dpid} - {desc} ---")
        for site in SITES:
            yrs = sorted({int(ym[:4]) for ym in target_months(dpid, site)})
            for yr in yrs:
                marker = os.path.join(SAVEPATH, f".done_{key}_{site}_{yr}")
                if os.path.exists(marker):
                    L(f"[skip] {key}/{site}/{yr} (already done)")
                    continue
                L(f"[get ] {key}/{site}/{yr} started {time.strftime('%H:%M:%S')}")
                try:
                    nu.by_file_aop(
                        dpid=dpid, site=site, year=str(yr),
                        include_provisional=True,
                        check_size=False,
                        savepath=SAVEPATH, token=TOKEN,
                        skip_if_exists=True,
                        overwrite="yes",
                    )
                    open(marker, "w").close()
                    L(f"[done] {key}/{site}/{yr}")
                except Exception as e:
                    L(f"[ERR ] {key}/{site}/{yr}: {type(e).__name__}: {e}")
    log.close()
    L(f"\nAll downloads completed: {time.strftime('%Y-%m-%d %H:%M:%S')}")


# ── Main ─────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["estimate", "download"], default="estimate", nargs="?")
    parser.add_argument("--product", type=str, default=None,
                        help="Filter by product key: VI, LAI, fPAR, ALB, CWI, CHM")
    args = parser.parse_args()

    if args.mode == "estimate":
        estimate(args.product)
    else:
        download(args.product)
