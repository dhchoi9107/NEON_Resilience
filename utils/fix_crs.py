"""
Fix missing or incorrect CRS on existing FSD GeoTIFF outputs.

Reads each *_FSD_10m.tif in the output directory, determines the correct
UTM zone from the NEON site code, and writes the CRS into the file header.

EPSG codes verified against NEON Data Portal API (data.neonscience.org/api/v0/sites)
and cross-checked with CRS extracted from LAZ file headers.

Usage:
  python fix_crs.py                # fix missing + correct wrong CRS
  python fix_crs.py --dry-run      # preview only, no changes
"""

import argparse
from pathlib import Path

import rasterio
from rasterio.crs import CRS

OUTPUT_DIR = Path("E:/neon_lidar/structural_diversity")

# NEON site -> EPSG (WGS 84 / UTM zone)
# Source: NEON Data Portal API (data.neonscience.org/api/v0/sites/{SITE})
# UTM zone = floor((longitude + 180) / 6) + 1, all Northern Hemisphere
SITE_EPSG = {
    # Northeast / Mid-Atlantic
    "HARV": 32618,  # UTM 18N  Harvard Forest, MA          lon=-72.173
    "BART": 32619,  # UTM 19N  Bartlett, NH                lon=-71.287
    "SCBI": 32617,  # UTM 17N  Smithsonian Conservation, VA lon=-78.139
    "SERC": 32618,  # UTM 18N  Smithsonian Env Research, MD lon=-76.560
    "BLAN": 32617,  # UTM 17N  Blandy Farm, VA             lon=-78.042
    # Southeast
    "ORNL": 32616,  # UTM 16N  Oak Ridge, TN               lon=-84.283
    "GRSM": 32617,  # UTM 17N  Great Smoky Mtns, TN/NC     lon=-83.502
    "MLBS": 32617,  # UTM 17N  Mountain Lake, VA           lon=-80.525
    "JERC": 32616,  # UTM 16N  Jones Ecological, GA        lon=-84.469
    "OSBS": 32617,  # UTM 17N  Ordway-Swisher, FL          lon=-81.993
    "TALL": 32616,  # UTM 16N  Talladega, AL               lon=-87.393
    "DELA": 32616,  # UTM 16N  Dead Lake, FL               lon=-87.804
    "DSNY": 32617,  # UTM 17N  Disney Wilderness, FL       lon=-81.436
    "LENO": 32616,  # UTM 16N  Lenoir Landing, AL          lon=-88.161
    # Great Lakes / Upper Midwest
    "UNDE": 32616,  # UTM 16N  UNDERC, MI                  lon=-89.537
    "STEI": 32616,  # UTM 16N  Steigerwaldt, WI            lon=-89.586
    "TREE": 32616,  # UTM 16N  Treehaven, WI               lon=-89.586
    "CHEQ": 32615,  # UTM 15N  Chequamegon-Nicolet, WI     (from LAZ headers)
    # Pacific Northwest
    "WREF": 32610,  # UTM 10N  Wind River, WA              lon=-121.952
    "ABBY": 32610,  # UTM 10N  Abby Road, WA               lon=-122.330
    # California / Sierra Nevada
    "SOAP": 32611,  # UTM 11N  Soaproot Saddle, CA         lon=-119.262
    "TEAK": 32611,  # UTM 11N  Lower Teakettle, CA         lon=-119.006
    "SJER": 32611,  # UTM 11N  San Joaquin, CA             lon=-119.732
    # Rocky Mountains
    "RMNP": 32613,  # UTM 13N  Rocky Mountain NP, CO       lon=-105.546
    "YELL": 32612,  # UTM 12N  Yellowstone, WY             lon=-110.539
    # Great Plains
    "KONZ": 32614,  # UTM 14N  Konza Prairie, KS           lon=-96.563
    "KONA": 32614,  # UTM 14N  Konza Ag, KS                lon=-96.613
    "WOOD": 32614,  # UTM 14N  Woodworth, ND               lon=-99.241
    "NOGP": 32614,  # UTM 14N  Northern Great Plains, ND   lon=-100.915
    "DCFS": 32614,  # UTM 14N  Dakota Coteau, ND           lon=-99.107
    "CPER": 32613,  # UTM 13N  Central Plains, CO          lon=-104.746
    "OAES": 32614,  # UTM 14N  Klemme Range, OK            lon=-99.059
    "STER": 32613,  # UTM 13N  North Sterling, CO          lon=-103.029
    # Alaska
    "BONA": 32606,  # UTM  6N  Caribou-Poker, AK           lon=-147.503
    "DEJU": 32606,  # UTM  6N  Delta Junction, AK          lon=-145.751
    "HEAL": 32606,  # UTM  6N  Healy, AK                   lon=-149.213
    # South-central
    "CLBJ": 32614,  # UTM 14N  LBJ Grassland, TX           lon=-97.570
    "UKFS": 32615,  # UTM 15N  UK Forest, KY               lon=-95.192
    # Puerto Rico
    "GUAN": 32619,  # UTM 19N  Guanica, PR                 lon=-66.869
    "LAJA": 32619,  # UTM 19N  Lajas, PR                   lon=-67.077
    # Hawaii
    "PUUM": 32605,  # UTM  5N  Pu'u Maka'ala, HI           lon=-155.317
}


def fix_file(tif_path: Path, dry_run: bool = False) -> str:
    """Fix CRS for a single GeoTIFF. Returns status string."""
    with rasterio.open(str(tif_path)) as src:
        existing_crs = src.crs

    # Extract site code from filename: {YEAR}_{SITE}_{N}_FSD_10m.tif
    stem = tif_path.stem
    parts = stem.split("_")
    site = parts[1] if len(parts) >= 2 else parts[0]

    if site not in SITE_EPSG:
        return f"  FAIL  {tif_path.name}  (unknown site: {site})"

    expected_epsg = SITE_EPSG[site]

    # Check existing CRS
    if existing_crs is not None and existing_crs.to_epsg() is not None:
        current_epsg = existing_crs.to_epsg()
        if current_epsg == expected_epsg:
            return f"  OK    {tif_path.name}  EPSG:{current_epsg}"
        # Wrong CRS - needs correction
        if dry_run:
            return f"  WOULD CORRECT  {tif_path.name}  EPSG:{current_epsg} -> {expected_epsg}"
        with rasterio.open(str(tif_path), "r+") as dst:
            dst.crs = CRS.from_epsg(expected_epsg)
        return f"  CORRECTED  {tif_path.name}  EPSG:{current_epsg} -> {expected_epsg}"

    # Missing CRS
    if dry_run:
        return f"  WOULD FIX  {tif_path.name}  -> EPSG:{expected_epsg}"

    with rasterio.open(str(tif_path), "r+") as dst:
        dst.crs = CRS.from_epsg(expected_epsg)

    return f"  FIXED  {tif_path.name}  -> EPSG:{expected_epsg}"


def main():
    parser = argparse.ArgumentParser(description="Fix missing/incorrect CRS on FSD GeoTIFFs")
    parser.add_argument("--dry-run", action="store_true", help="Preview only, no changes")
    args = parser.parse_args()

    tifs = sorted(OUTPUT_DIR.glob("*_FSD_10m.tif"))
    if not tifs:
        print(f"No *_FSD_10m.tif files found in {OUTPUT_DIR}")
        return

    print(f"Found {len(tifs)} FSD files in {OUTPUT_DIR}")
    if args.dry_run:
        print("(DRY RUN - no files will be modified)\n")

    ok = fixed = corrected = failed = 0
    for tif in tifs:
        status = fix_file(tif, dry_run=args.dry_run)
        print(status)
        if "OK" in status:
            ok += 1
        elif "FIXED" in status or "WOULD FIX" in status:
            fixed += 1
        elif "CORRECTED" in status or "WOULD CORRECT" in status:
            corrected += 1
        else:
            failed += 1

    print(f"\nDone: {ok} ok, {fixed} fixed, {corrected} corrected, {failed} failed")


if __name__ == "__main__":
    main()
