"""
NEON Vegetation Structure (DP1.10098.001) Download
===================================================
Downloads tree species inventory data (CSV tables) for all curated sites.
This is observational (IS) data, not AOP — uses neonutilities.load_by_product().

Key tables:
  vst_mappingandtagging   - species identity, stem location (UTM)
  vst_apparentindividual  - annual measurements (diameter, height, status)
  vst_perplotperyear      - plot-level metadata

Usage:
  python neon_veg_structure_download.py estimate    # check availability
  python neon_veg_structure_download.py download     # download all
  python neon_veg_structure_download.py download --site HARV  # single site
"""

import os
import sys
import argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from site_config import (
    SITES, YEARS, VEG_STRUCT_DIR, TOKEN, NEON_API, api_get, available_months,
)

DPID = "DP1.10098.001"


def estimate():
    """Check data availability per site."""
    print(f"Product: {DPID} (Vegetation Structure)")
    print(f"Sites: {len(SITES)}")
    print(f"Output: {VEG_STRUCT_DIR}\n")

    total_months = 0
    for site in SITES:
        months = available_months(DPID, site)
        year_months = [m for m in months if any(m.startswith(str(y)) for y in YEARS)]
        n = len(year_months)
        total_months += n
        yr_range = ""
        if year_months:
            yr_range = f"{year_months[0][:4]}-{year_months[-1][:4]}"
        print(f"  {site}: {n:3d} months  ({yr_range})")

    print(f"\nTotal: {total_months} site-months across {len(SITES)} sites")


def download(site_filter=None):
    """Download vegetation structure data using neonutilities."""
    try:
        import neonutilities
    except ImportError:
        print("ERROR: neonutilities not installed.")
        print("  pip install neonutilities")
        sys.exit(1)

    VEG_STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    sites = [site_filter] if site_filter else SITES
    start_date = f"{min(YEARS)}-01"
    end_date = f"{max(YEARS)}-12"

    for site in sites:
        marker = VEG_STRUCT_DIR / f".done_{site}"
        if marker.exists():
            print(f"  SKIP {site} (already downloaded)")
            continue

        print(f"\n{'='*60}")
        print(f"  Downloading {DPID} for {site} ({start_date} to {end_date})")
        print(f"{'='*60}")

        try:
            result = neonutilities.load_by_product(
                dpid=DPID,
                site=site,
                startdate=start_date,
                enddate=end_date,
                include_provisional=True,
                token=TOKEN or "",
                check_size=False,
            )
            # result is a dict of {table_name: DataFrame}
            # Save each table as CSV
            import pandas as pd
            for table_name, df in result.items():
                if isinstance(df, pd.DataFrame) and not df.empty:
                    out_csv = VEG_STRUCT_DIR / f"{table_name}.csv"
                    if out_csv.exists():
                        existing = pd.read_csv(out_csv, low_memory=False)
                        df = pd.concat([existing, df], ignore_index=True).drop_duplicates()
                    df.to_csv(out_csv, index=False)
            marker.touch()
            print(f"  DONE {site} ({len(result)} tables)")

        except Exception as e:
            print(f"  ERROR {site}: {e}")
            continue


def main():
    parser = argparse.ArgumentParser(description="Download NEON vegetation structure data")
    parser.add_argument("mode", choices=["estimate", "download"],
                        help="estimate: check availability; download: fetch data")
    parser.add_argument("--site", type=str, default=None,
                        help="Download only this site (e.g., HARV)")
    args = parser.parse_args()

    if args.mode == "estimate":
        estimate()
    else:
        download(site_filter=args.site)


if __name__ == "__main__":
    main()
