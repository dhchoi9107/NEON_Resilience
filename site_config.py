"""
Shared configuration for NEON Functional Biodiversity analysis pipeline.
========================================================================
Centralizes site lists, CRS mappings, NEON domain codes, API helpers,
and path constants used by all downstream scripts.
"""

import os
import json
import time
import urllib.request
import urllib.error
from pathlib import Path

# ─── Paths ──────────────────────────────────────────────────────────────────

NEON_BASE   = Path("E:/neon_lidar/DP1.30003.001")
FSD_DIR     = Path("E:/neon_lidar/structural_diversity")
VI_DIR      = Path("E:/neon_lidar/vegetation_indices")
VEG_STRUCT_DIR = Path("E:/neon_lidar/vegetation_structure")
TAX_DIV_DIR = Path("E:/neon_lidar/taxonomic_diversity")
SPEC_DIV_DIR = Path("E:/neon_lidar/spectral_diversity")
FUNC_DIV_DIR = Path("E:/neon_lidar/functional_diversity")
ENV_HET_DIR = Path("E:/neon_lidar/env_heterogeneity")
DHI_DIR     = Path("E:/neon_lidar/productivity_dhi")
MODEL_DIR   = Path("E:/neon_lidar/model_results")

CELL_SIZE = 10  # Sentinel-2 aligned grid resolution (meters)

# ─── Sites (curated temperate forest set) ───────────────────────────────────

SITES = sorted([
    "HARV", "BART", "UNDE", "BLAN", "SCBI", "SERC", "MLBS", "ORNL", "GRSM",
    "WREF", "ABBY", "SOAP", "TEAK", "RMNP", "JERC", "OSBS", "TALL", "STEI",
    "TREE",
])

YEARS = list(range(2013, 2027))

# Override via environment
if os.environ.get("SITES_ENV"):
    SITES = [s.strip() for s in os.environ["SITES_ENV"].split(",") if s.strip()]
if os.environ.get("YEARS_ENV"):
    YEARS = [int(y.strip()) for y in os.environ["YEARS_ENV"].split(",") if y.strip()]

# ─── NEON Domain codes (for mixed-effects model random effects) ─────────────
# Source: NEON Data Portal API (data.neonscience.org/api/v0/sites/{SITE})

SITE_DOMAIN = {
    "HARV": "D01",  # Northeast
    "BART": "D01",  # Northeast
    "SCBI": "D02",  # Mid-Atlantic
    "SERC": "D02",  # Mid-Atlantic
    "BLAN": "D02",  # Mid-Atlantic
    "ORNL": "D07",  # Appalachians & Cumberland Plateau
    "GRSM": "D07",  # Appalachians & Cumberland Plateau
    "MLBS": "D07",  # Appalachians & Cumberland Plateau
    "JERC": "D03",  # Southeast
    "OSBS": "D03",  # Southeast
    "TALL": "D08",  # Ozarks Complex
    "UNDE": "D05",  # Great Lakes
    "STEI": "D05",  # Great Lakes
    "TREE": "D05",  # Great Lakes
    "CHEQ": "D05",  # Great Lakes
    "WREF": "D16",  # Pacific Northwest
    "ABBY": "D16",  # Pacific Northwest
    "SOAP": "D17",  # Pacific Southwest
    "TEAK": "D17",  # Pacific Southwest
    "RMNP": "D10",  # Central Plains & Northern Rockies
}

# ─── Site EPSG codes (WGS 84 / UTM zone) ───────────────────────────────────
# Source: NEON Data Portal API coordinates, verified against LAZ headers

SITE_EPSG = {
    "HARV": 32618,  # UTM 18N  lon=-72.173
    "BART": 32619,  # UTM 19N  lon=-71.287
    "SCBI": 32617,  # UTM 17N  lon=-78.139
    "SERC": 32618,  # UTM 18N  lon=-76.560
    "BLAN": 32617,  # UTM 17N  lon=-78.042
    "ORNL": 32616,  # UTM 16N  lon=-84.283
    "GRSM": 32617,  # UTM 17N  lon=-83.502
    "MLBS": 32617,  # UTM 17N  lon=-80.525
    "JERC": 32616,  # UTM 16N  lon=-84.469
    "OSBS": 32617,  # UTM 17N  lon=-81.993
    "TALL": 32616,  # UTM 16N  lon=-87.393
    "UNDE": 32616,  # UTM 16N  lon=-89.537
    "STEI": 32616,  # UTM 16N  lon=-89.586
    "TREE": 32616,  # UTM 16N  lon=-89.586
    "CHEQ": 32615,  # UTM 15N  (from LAZ headers)
    "WREF": 32610,  # UTM 10N  lon=-121.952
    "ABBY": 32610,  # UTM 10N  lon=-122.330
    "SOAP": 32611,  # UTM 11N  lon=-119.262
    "TEAK": 32611,  # UTM 11N  lon=-119.006
    "RMNP": 32613,  # UTM 13N  lon=-105.546
}

# ─── Site coordinates (for plot-level spatial referencing) ──────────────────
# Source: NEON Data Portal API

SITE_COORDS = {
    "HARV": (42.5369, -72.1727),
    "BART": (44.0639, -71.2874),
    "SCBI": (38.8929, -78.1395),
    "SERC": (38.8901, -76.5600),
    "BLAN": (39.0337, -78.0418),
    "ORNL": (35.9641, -84.2826),
    "GRSM": (35.6890, -83.5020),
    "MLBS": (37.3783, -80.5248),
    "JERC": (31.1948, -84.4686),
    "OSBS": (29.6893, -81.9934),
    "TALL": (32.9505, -87.3933),
    "UNDE": (46.2339, -89.5373),
    "STEI": (45.5089, -89.5864),
    "TREE": (45.4937, -89.5857),
    "WREF": (45.8205, -121.9519),
    "ABBY": (45.7624, -122.3303),
    "SOAP": (37.0334, -119.2622),
    "TEAK": (37.0058, -119.0060),
    "RMNP": (40.2759, -105.5460),
}

# ─── NEON API ───────────────────────────────────────────────────────────────

NEON_API = "https://data.neonscience.org/api/v0"


def load_token():
    """Load NEON API token from env var or local files."""
    token = os.environ.get("NEON_TOKEN")
    if not token:
        for tf in (
            os.path.join(os.path.dirname(os.path.abspath(__file__)), ".neon_token"),
            os.path.expanduser("~/.neon_token"),
            "E:/neon_lidar/_code/.neon_token",
        ):
            if os.path.exists(tf):
                token = open(tf).read().strip()
                break
    return token


TOKEN = load_token()


def api_get(url, token=None):
    """GET a NEON API endpoint with retry and optional auth token."""
    t = token or TOKEN
    req = urllib.request.Request(url)
    if t:
        req.add_header("X-API-Token", t)
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
    """Return list of available months for a product at a site."""
    if dpid not in _PROD_CACHE:
        d = api_get(f"{NEON_API}/products/{dpid}")["data"]
        _PROD_CACHE[dpid] = {}
        for s in d["siteCodes"]:
            _PROD_CACHE[dpid][s["siteCode"]] = sorted(s.get("availableMonths", []))
    return _PROD_CACHE[dpid].get(site, [])


# ─── FSD file discovery ─────────────────────────────────────────────────────

def get_fsd_files(site=None):
    """Return sorted list of FSD GeoTIFF paths, optionally filtered by site.

    File naming convention: {YEAR}_{SITE}_{N}_FSD_10m.tif
    Returns list of (year, site, path) tuples sorted by year.
    """
    results = []
    for f in FSD_DIR.glob("*_FSD_10m.tif"):
        parts = f.stem.split("_")
        if len(parts) < 3:
            continue
        yr = int(parts[0])
        st = parts[1]
        if site and st != site:
            continue
        results.append((yr, st, f))
    return sorted(results)


def get_fsd_band_index(band_name):
    """Return 1-based band index for a given FSD metric name."""
    return FSD_BANDS.index(band_name) + 1


# FSD band names (must match compute_fsd.py BAND_NAMES)
FSD_BANDS = [
    "rumple", "top_rugosity", "mean_max_canopy_ht", "max_canopy_ht",
    "deepgap_fraction", "meanH", "vert_sd", "vertCV", "mean_sd", "sd_sd",
    "GC", "GFP", "VCI", "q25", "q50", "q75", "q95", "HeightRatio",
    "FHD", "LAI", "LAI_subcanopy",
]


def ensure_dirs():
    """Create all output directories if they don't exist."""
    for d in [TAX_DIV_DIR, SPEC_DIV_DIR, FUNC_DIV_DIR,
              ENV_HET_DIR, DHI_DIR, MODEL_DIR, VEG_STRUCT_DIR]:
        d.mkdir(parents=True, exist_ok=True)
