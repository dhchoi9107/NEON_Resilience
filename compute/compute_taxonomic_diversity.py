"""
Taxonomic Diversity Calculator from NEON Vegetation Structure
==============================================================
Computes alpha, beta, and gamma taxonomic diversity from NEON
vegetation structure data (DP1.10098.001).

Alpha (per plot, per year):
  - Species richness (S)
  - Shannon diversity (H')
  - Simpson diversity (1-D)
  - Stem abundance (N)

Beta (pairwise between plots within site, per year):
  - Bray-Curtis dissimilarity
  - Jaccard dissimilarity
  - Sorensen dissimilarity
  - Baselga (2010) partitioning: turnover (beta_sim) + nestedness (beta_sne)

Gamma (per site, per year):
  - Total species richness across all plots

Usage:
  python compute_taxonomic_diversity.py                  # all sites
  python compute_taxonomic_diversity.py --site HARV      # single site
"""

import argparse
import sys
import warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis, pdist, squareform

from site_config import SITES, VEG_STRUCT_DIR, TAX_DIV_DIR

warnings.filterwarnings("ignore")


# ─── Data Loading ────────────────────────────────────────────────────────────

def find_vst_tables(base_dir):
    """Find stacked CSV tables from neonutilities download."""
    base = Path(base_dir)

    # neonutilities creates: filesToStack/stackedFiles/ with stacked CSVs
    # Or directly in the download folder depending on version
    candidates = [
        base / "filesToStack" / "stackedFiles",
        base,
    ]

    # Also search recursively for the stacked files
    for pattern in ["**/vst_apparentindividual.csv", "**/stackedFiles/vst_apparentindividual.csv"]:
        found = list(base.glob(pattern))
        if found:
            candidates.insert(0, found[0].parent)

    tables = {}
    for d in candidates:
        if not d.exists():
            continue
        for name in ["vst_mappingandtagging", "vst_apparentindividual", "vst_perplotperyear"]:
            path = d / f"{name}.csv"
            if path.exists() and name not in tables:
                tables[name] = path

    return tables


def load_vst_data(base_dir):
    """Load and join vegetation structure tables into a clean dataframe.

    Returns DataFrame with columns:
      siteID, plotID, year, taxonID, stemDiameter, height, plantStatus,
      pointID, stemTag, easting, northing
    """
    tables = find_vst_tables(base_dir)
    required = ["vst_mappingandtagging", "vst_apparentindividual"]
    for t in required:
        if t not in tables:
            raise FileNotFoundError(
                f"Cannot find {t}.csv in {base_dir}. "
                f"Run neon_veg_structure_download.py first."
            )

    print(f"  Loading: {tables['vst_mappingandtagging']}")
    tag = pd.read_csv(tables["vst_mappingandtagging"], low_memory=False, usecols=[
        "individualID", "siteID", "plotID", "taxonID",
        "pointID",
    ]).drop_duplicates(subset=["individualID"])

    print(f"  Loading: {tables['vst_apparentindividual']}")
    ind = pd.read_csv(tables["vst_apparentindividual"], low_memory=False, usecols=[
        "individualID", "siteID", "plotID", "date",
        "stemDiameter", "height", "plantStatus",
    ])

    # Extract year from date
    ind["date"] = pd.to_datetime(ind["date"], errors="coerce")
    ind["year"] = ind["date"].dt.year.astype("Int64")
    ind = ind.dropna(subset=["year"])
    ind["year"] = ind["year"].astype(int)

    # Join species info
    df = ind.merge(tag[["individualID", "taxonID", "pointID"]],
                   on="individualID", how="left")

    # Filter to live stems only
    df = df[df["plantStatus"].str.contains("Live", case=False, na=False)]

    # Filter to woody stems with DBH measurement (exclude shrubs/seedlings)
    n_before = len(df)
    df = df[df["stemDiameter"].notna() & (df["stemDiameter"] >= 10.0)]
    print(f"  DBH filter (>=10cm): {n_before:,} -> {len(df):,} "
          f"({len(df)/n_before*100:.0f}% retained)")

    # Remove unknown/unidentified taxa
    df = df.dropna(subset=["taxonID"])
    df = df[~df["taxonID"].str.contains("Unknown|unidentified|spp\\.", case=False, na=False)]

    # Load plot coordinates if available
    if "vst_perplotperyear" in tables:
        ppy = pd.read_csv(tables["vst_perplotperyear"], low_memory=False, usecols=[
            "siteID", "plotID", "totalSampledAreaTrees", "nlcdClass",
        ]).drop_duplicates(subset=["siteID", "plotID"])
        df = df.merge(ppy, on=["siteID", "plotID"], how="left")

    print(f"  Records: {len(df):,} live stems, "
          f"{df['taxonID'].nunique()} species, "
          f"{df['plotID'].nunique()} plots")

    return df


# ─── Diversity Metrics ───────────────────────────────────────────────────────

def shannon(counts):
    """Shannon diversity index H' = -sum(p_i * ln(p_i))."""
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0.0
    p = counts / counts.sum()
    return -np.sum(p * np.log(p))


def simpson(counts):
    """Simpson diversity index 1-D = 1 - sum(p_i^2)."""
    counts = counts[counts > 0]
    if len(counts) <= 1:
        return 0.0
    p = counts / counts.sum()
    return 1.0 - np.sum(p ** 2)


def compute_alpha(df, site, year):
    """Compute alpha diversity metrics per plot."""
    sub = df[(df["siteID"] == site) & (df["year"] == year)]
    if sub.empty:
        return pd.DataFrame()

    results = []
    for plot_id, grp in sub.groupby("plotID"):
        counts = grp.groupby("taxonID").size().values
        results.append({
            "siteID": site,
            "plotID": plot_id,
            "year": year,
            "richness": len(counts),
            "abundance": int(counts.sum()),
            "shannon": shannon(counts),
            "simpson": simpson(counts),
        })

    return pd.DataFrame(results)


def baselga_partition(pa_matrix):
    """Baselga (2010) partitioning of Sorensen beta diversity.

    pa_matrix: presence-absence matrix (plots x species), binary.
    Returns dict with site-level mean beta_sor, beta_sim (turnover),
    beta_sne (nestedness).
    """
    n_plots = pa_matrix.shape[0]
    if n_plots < 2:
        return {"beta_sor": np.nan, "beta_sim": np.nan, "beta_sne": np.nan}

    sor_vals = []
    sim_vals = []
    sne_vals = []

    for i in range(n_plots):
        for j in range(i + 1, n_plots):
            a = np.sum((pa_matrix[i] > 0) & (pa_matrix[j] > 0))  # shared
            b = np.sum((pa_matrix[i] > 0) & (pa_matrix[j] == 0))  # only in i
            c = np.sum((pa_matrix[i] == 0) & (pa_matrix[j] > 0))  # only in j

            denom_sor = 2 * a + b + c
            if denom_sor == 0:
                continue

            beta_sor = (b + c) / denom_sor
            denom_sim = a + min(b, c)
            beta_sim = min(b, c) / denom_sim if denom_sim > 0 else 0.0
            beta_sne = beta_sor - beta_sim

            sor_vals.append(beta_sor)
            sim_vals.append(beta_sim)
            sne_vals.append(beta_sne)

    return {
        "beta_sor": np.mean(sor_vals) if sor_vals else np.nan,
        "beta_sim": np.mean(sim_vals) if sim_vals else np.nan,
        "beta_sne": np.mean(sne_vals) if sne_vals else np.nan,
    }


def compute_beta(df, site, year):
    """Compute beta diversity metrics between plots at a site-year."""
    sub = df[(df["siteID"] == site) & (df["year"] == year)]
    if sub.empty:
        return {}

    # Build abundance matrix (plots x species)
    pivot = sub.groupby(["plotID", "taxonID"]).size().unstack(fill_value=0)
    if pivot.shape[0] < 2:
        return {"siteID": site, "year": year, "n_plots": pivot.shape[0],
                "bray_mean": np.nan, "jaccard_mean": np.nan,
                "beta_sor": np.nan, "beta_sim": np.nan, "beta_sne": np.nan}

    abundance = pivot.values.astype(float)
    pa = (abundance > 0).astype(float)

    # Bray-Curtis (abundance-based)
    bray_dists = pdist(abundance, metric="braycurtis")
    bray_mean = np.mean(bray_dists)

    # Jaccard (presence-absence)
    jaccard_dists = pdist(pa, metric="jaccard")
    jaccard_mean = np.mean(jaccard_dists)

    # Baselga partitioning
    baselga = baselga_partition(pa)

    return {
        "siteID": site,
        "year": year,
        "n_plots": pivot.shape[0],
        "n_species": pivot.shape[1],
        "bray_mean": float(bray_mean),
        "jaccard_mean": float(jaccard_mean),
        **baselga,
    }


def compute_gamma(df, site, year):
    """Compute gamma diversity (site-level total richness)."""
    sub = df[(df["siteID"] == site) & (df["year"] == year)]
    if sub.empty:
        return {}
    return {
        "siteID": site,
        "year": year,
        "gamma_richness": sub["taxonID"].nunique(),
        "gamma_abundance": len(sub),
    }


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compute taxonomic diversity from NEON VST data")
    parser.add_argument("--site", type=str, default=None, help="Process single site")
    args = parser.parse_args()

    TAX_DIV_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading NEON vegetation structure data...")
    df = load_vst_data(VEG_STRUCT_DIR)

    available_sites = sorted(df["siteID"].unique())
    target_sites = [args.site] if args.site else [s for s in SITES if s in available_sites]

    if not target_sites:
        print("No sites found in the data. Available:", available_sites)
        return

    print(f"\nProcessing {len(target_sites)} sites: {target_sites}")

    # ── Per-year diversity (kept for temporal analysis) ──
    all_alpha = []
    all_beta = []
    all_gamma = []

    for site in target_sites:
        site_df = df[df["siteID"] == site]
        site_years = sorted(site_df["year"].unique())
        print(f"\n{'='*50}")
        print(f"  {site}: {len(site_years)} years ({min(site_years)}-{max(site_years)})")

        for year in site_years:
            alpha = compute_alpha(df, site, year)
            if not alpha.empty:
                all_alpha.append(alpha)
            beta = compute_beta(df, site, year)
            if beta:
                all_beta.append(beta)
            gamma = compute_gamma(df, site, year)
            if gamma:
                all_gamma.append(gamma)

        n_plots = site_df["plotID"].nunique()
        n_species = site_df["taxonID"].nunique()
        print(f"    {n_plots} plots, {n_species} species")

    # ── Pooled diversity (all years combined per plot) ──
    # Trees are long-lived — pooling across years captures the true species
    # composition better than any single year's partial survey.
    print(f"\n{'='*50}")
    print("Computing POOLED diversity (all years combined)...")

    pooled_alpha = []
    pooled_beta = []
    pooled_gamma = []

    for site in target_sites:
        site_df = df[df["siteID"] == site]
        # Deduplicate: one record per individual per plot (keep unique individuals)
        site_pooled = site_df.drop_duplicates(subset=["plotID", "individualID"])

        # Alpha per plot (pooled)
        for plot_id, grp in site_pooled.groupby("plotID"):
            counts = grp.groupby("taxonID").size().values
            pooled_alpha.append({
                "siteID": site,
                "plotID": plot_id,
                "year": "pooled",
                "richness": len(counts),
                "abundance": int(counts.sum()),
                "shannon": shannon(counts),
                "simpson": simpson(counts),
            })

        # Beta (pooled)
        pivot = site_pooled.groupby(["plotID", "taxonID"]).size().unstack(fill_value=0)
        if pivot.shape[0] >= 2:
            abundance = pivot.values.astype(float)
            pa = (abundance > 0).astype(float)
            from scipy.spatial.distance import pdist
            bray_dists = pdist(abundance, metric="braycurtis")
            jaccard_dists = pdist(pa, metric="jaccard")
            baselga = baselga_partition(pa)
            pooled_beta.append({
                "siteID": site,
                "year": "pooled",
                "n_plots": pivot.shape[0],
                "n_species": pivot.shape[1],
                "bray_mean": float(np.mean(bray_dists)),
                "jaccard_mean": float(np.mean(jaccard_dists)),
                **baselga,
            })

        # Gamma (pooled)
        pooled_gamma.append({
            "siteID": site,
            "year": "pooled",
            "gamma_richness": site_pooled["taxonID"].nunique(),
            "gamma_abundance": len(site_pooled),
        })

        print(f"  {site}: {pivot.shape[0] if pivot.shape[0] >= 2 else 0} plots, "
              f"{site_pooled['taxonID'].nunique()} species (pooled)")

    # ── Save outputs ──
    # Per-year
    if all_alpha:
        alpha_df = pd.concat(all_alpha, ignore_index=True)
        alpha_df.to_csv(TAX_DIV_DIR / "alpha_diversity.csv", index=False)
        print(f"\nAlpha (per-year): {len(alpha_df)} plot-years")

    if all_beta:
        beta_df = pd.DataFrame(all_beta)
        beta_df.to_csv(TAX_DIV_DIR / "beta_diversity.csv", index=False)
        print(f"Beta (per-year): {len(beta_df)} site-years")

    if all_gamma:
        gamma_df = pd.DataFrame(all_gamma)
        gamma_df.to_csv(TAX_DIV_DIR / "gamma_diversity.csv", index=False)
        print(f"Gamma (per-year): {len(gamma_df)} site-years")

    # Pooled
    if pooled_alpha:
        palpha_df = pd.DataFrame(pooled_alpha)
        palpha_df.to_csv(TAX_DIV_DIR / "alpha_diversity_pooled.csv", index=False)
        print(f"Alpha (pooled): {len(palpha_df)} plots")

    if pooled_beta:
        pbeta_df = pd.DataFrame(pooled_beta)
        pbeta_df.to_csv(TAX_DIV_DIR / "beta_diversity_pooled.csv", index=False)
        print(f"Beta (pooled): {len(pbeta_df)} sites")

    if pooled_gamma:
        pgamma_df = pd.DataFrame(pooled_gamma)
        pgamma_df.to_csv(TAX_DIV_DIR / "gamma_diversity_pooled.csv", index=False)
        print(f"Gamma (pooled): {len(pgamma_df)} sites")

    # Combined for modeling (pooled version — primary)
    if pooled_alpha and pooled_beta:
        palpha_df = pd.DataFrame(pooled_alpha)
        pbeta_df = pd.DataFrame(pooled_beta)
        pgamma_df = pd.DataFrame(pooled_gamma)

        alpha_summary = palpha_df.groupby("siteID").agg({
            "richness": "mean",
            "abundance": "mean",
            "shannon": "mean",
            "simpson": "mean",
        }).reset_index()
        alpha_summary.columns = ["siteID",
                                 "alpha_richness_mean", "alpha_abundance_mean",
                                 "alpha_shannon_mean", "alpha_simpson_mean"]

        combined = pbeta_df.drop(columns=["year"]).merge(alpha_summary, on="siteID", how="outer")
        combined = combined.merge(pgamma_df.drop(columns=["year"]), on="siteID", how="outer")

        combined.to_csv(TAX_DIV_DIR / "taxonomic_pooled.csv", index=False)
        print(f"Combined (pooled): {TAX_DIV_DIR / 'taxonomic_pooled.csv'}")


if __name__ == "__main__":
    main()
