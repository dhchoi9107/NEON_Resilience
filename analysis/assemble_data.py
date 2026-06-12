"""
Assemble All Diversity & Productivity Data for Modeling
========================================================
Merges outputs from all pipeline components into a single analysis dataset.

Inputs:
  - taxonomic_diversity/  (alpha, beta, gamma CSVs)
  - spectral_diversity/   (spectral_diversity_all.csv)
  - functional_diversity/ (functional_diversity_all.csv)
  - env_heterogeneity/    (heterogeneity_all.csv)
  - productivity_dhi/     (productivity_summary.csv)

Output:
  - model_results/assembled_data.csv  (ready for R lme4)

Usage:
  python analysis/assemble_data.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
import numpy as np

from site_config import (
    TAX_DIV_DIR, SPEC_DIV_DIR, FUNC_DIV_DIR, ENV_HET_DIR, DHI_DIR,
    MODEL_DIR, SITE_DOMAIN, SITE_COORDS,
)


def load_csv(path, label):
    """Load CSV if exists, else return empty DataFrame."""
    if path.exists():
        df = pd.read_csv(path)
        print(f"  {label}: {len(df)} rows, {list(df.columns[:5])}...")
        return df
    else:
        print(f"  {label}: NOT FOUND ({path})")
        return pd.DataFrame()


def main():
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    print("Assembling analysis dataset\n")

    # ── 1) Taxonomic diversity (pooled = primary) ──
    tax_pooled = load_csv(TAX_DIV_DIR / "taxonomic_pooled.csv", "Taxonomic (pooled)")
    tax_alpha_pooled = load_csv(TAX_DIV_DIR / "alpha_diversity_pooled.csv", "Alpha (pooled)")
    # Also keep per-year for temporal analysis
    tax_combined = load_csv(TAX_DIV_DIR / "taxonomic_all.csv", "Taxonomic (per-year)")
    tax_alpha = load_csv(TAX_DIV_DIR / "alpha_diversity.csv", "Alpha (per-year)")

    # ── 2) Spectral diversity ──
    spec = load_csv(SPEC_DIV_DIR / "spectral_diversity_all.csv", "Spectral diversity")

    # ── 3) Functional diversity ──
    func = load_csv(FUNC_DIV_DIR / "functional_diversity_all.csv", "Functional diversity")

    # ── 4) Environmental heterogeneity ──
    het = load_csv(ENV_HET_DIR / "heterogeneity_all.csv", "Heterogeneity")

    # ── 5) Productivity ──
    dhi = load_csv(DHI_DIR / "productivity_summary.csv", "Productivity (DHI)")
    if "site" in dhi.columns:
        dhi = dhi.rename(columns={"site": "siteID"})

    # ── Merge: site-year level ──
    print("\nMerging datasets...")

    # Start with taxonomic (site-year level summary)
    if not tax_combined.empty:
        base = tax_combined.copy()
    elif not tax_alpha.empty:
        # Aggregate alpha to site-year
        base = tax_alpha.groupby(["siteID", "year"]).agg({
            "richness": "mean",
            "shannon": "mean",
            "simpson": "mean",
            "abundance": "mean",
        }).reset_index()
    else:
        # Build base from heterogeneity or other available data
        frames = [df for df in [spec, het] if not df.empty and "siteID" in df.columns]
        if frames:
            base = frames[0][["siteID", "year"]].drop_duplicates()
        else:
            print("ERROR: No data available to build base dataset.")
            return

    # Merge spectral
    if not spec.empty and "siteID" in spec.columns:
        base = base.merge(spec, on=["siteID", "year"], how="outer", suffixes=("", "_spec"))

    # Merge heterogeneity
    if not het.empty and "siteID" in het.columns:
        base = base.merge(het, on=["siteID", "year"], how="left", suffixes=("", "_het"))

    # Merge DHI (site-level, no year)
    if not dhi.empty and "siteID" in dhi.columns:
        dhi_cols = [c for c in dhi.columns if c != "year"]
        base = base.merge(dhi[dhi_cols].drop_duplicates(subset=["siteID"]),
                          on="siteID", how="left", suffixes=("", "_dhi"))

    # Merge functional diversity (plot-level → aggregate to site-year)
    if not func.empty and "siteID" in func.columns:
        func_agg = func.groupby(["siteID", "year"]).agg({
            "func_FRic": "mean",
            "func_FDiv": "mean",
            "func_FEve": "mean",
            "func_RaoQ": "mean",
            "func_beta_rao": "first",
        }).reset_index()
        base = base.merge(func_agg, on=["siteID", "year"], how="left",
                          suffixes=("", "_func"))

    # ── Add metadata ──
    base["domain"] = base["siteID"].map(SITE_DOMAIN)
    base["latitude"] = base["siteID"].map(lambda s: SITE_COORDS.get(s, (np.nan,))[0])
    base["longitude"] = base["siteID"].map(lambda s: SITE_COORDS.get(s, (np.nan, np.nan))[1])

    # ── Report ──
    print(f"\nAssembled dataset: {len(base)} rows, {len(base.columns)} columns")
    print(f"Sites: {base['siteID'].nunique()}")
    if "year" in base.columns:
        print(f"Years: {base['year'].min()}-{base['year'].max()}")

    # Missing data report
    n_missing = base.isnull().sum()
    cols_with_missing = n_missing[n_missing > 0]
    if len(cols_with_missing) > 0:
        print(f"\nColumns with missing values ({len(cols_with_missing)}):")
        for col, n in cols_with_missing.items():
            pct = n / len(base) * 100
            print(f"  {col}: {n} ({pct:.0f}%)")

    # ── Save ──
    out_path = MODEL_DIR / "assembled_data.csv"
    base.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    # Also save plot-level alpha diversity separately for plot-level models
    if not tax_alpha.empty:
        alpha_path = MODEL_DIR / "plot_alpha_diversity.csv"
        tax_alpha.to_csv(alpha_path, index=False)
        print(f"Saved: {alpha_path}")

    # Save functional diversity (plot-level) for plot-level models
    if not func.empty:
        func_path = MODEL_DIR / "plot_functional_diversity.csv"
        func.to_csv(func_path, index=False)
        print(f"Saved: {func_path}")


if __name__ == "__main__":
    main()
