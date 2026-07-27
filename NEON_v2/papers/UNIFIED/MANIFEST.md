# UNIFIED paper — figure & result manifest

Single integrative paper. Files copied from `NEON_v2/figures/` and `NEON_v2/results/` (originals untouched). Organized by objective. **Artifact analyses excluded** (see bottom).

Title: *Static canopy state, structural dynamics, and ecosystem productivity encode distinct dimensions of forest biodiversity across North America.*

## O0 — Unified statistical framework (§2.5) — cross-cutting
**figures/O0_framework/O0_framework.png** ★ (script `154_framework_figure.py`): (A) incremental R² beyond domain State→+Dynamics→+Function with LRT stars; (B) unique/shared variance partition — structure dominant, spectral≈0 for beta.
**results/O0_framework/** (script `scripts/obj23/153_nested_framework.py`; see its README). VIF<5 screening → nested M0 domain→M1 state→M2 +dynamics→M3 +function (OLS+domain+site-clustered SE) with ΔR²/ΔAIC/LRT → unique/shared variance partition → interaction FDR → Moran's I + spatial-robust refit. Headline: dynamics add beyond state (all responses, p≤.001); GPP adds beyond state+dynamics (alpha p≤.03, turnover p=2e-5; nestedness n.s.); spectral unique≈0 for beta; residual autocorrelation present but inferences robust; 20/68 interactions survive BH.

---

## O1 — Static canopy state (state) → alpha/beta complementarity
**Finding**: spectral VI strongest for alpha (SAVI/EVI β≈+0.38); static structure uniquely carries beta (LCBD turnover). Structure & spectra complementary.
- **figures/O1_static/**: `K01` (headline RS→diversity, 4 responses) ★ · `F01*` forest · `F02*` scatter · `F03` feature effect · `F04` variance decomposition (domain>site>plot) · `F05*` by domain · `F06` by site · `Fig_v2_forest`
- **results/O1_static/**: `v2_coeff.csv` · `v2_variance_decomp.csv` · `obj1_specdiv.csv` (archived combined spectral-diversity — reference only, paper uses individual VIs)

## O2 — Structural dynamics (process) → beta increment + succession + disturbance
**Finding**: (a) temporal structural trends add beta info beyond mean (19/60 pairs, all structural). (b) age↑→structural change↓ (VCI/LAI trend β=−0.25/−0.22; old-growth Deep_Gap +0.28 = aggradation→gap-phase; partial-regression survivors = VCI/LAI/Deep_Gap/Canopy_Ht). (c) press (insect ΔLAI −0.18) vs pulse (fire/windthrow +0.22). (d) complexity buffering real only for Deep_Gap/Gini after RTM control.
- **figures/O2_dynamics/**: `M01` trend increment ★ · `F07` timeseries · `F08` trend distributions · `F09` coupling heatmaps · `N01`–`N08` age↔structural change (`N04` partial regression ★, `N06` Canopy_Ht Simpson) · `N11` 20-yr binned trends (robustness) · `J01` press/pulse ★ · `J02`–`J06` complexity buffer / RTM control
- **results/O2_dynamics/**: `trend_increment.csv` · `stand_age_trendage_{gami,30m}.csv` · `stand_age_models_{gami,30m}.csv` · `stand_age_crossval.csv` · `age_binned_trend.csv` · `choi2023.csv` (press/pulse) · `complexity_rtm.csv`

## O3 — Ecosystem productivity (function) → species–energy
**Finding**: independent GPP → diversity monotonic positive. MODIS β=+0.22, PML β=+0.17; tower-GPP correlations MODIS 0.67 / PML 0.69. Species–energy, not unimodal.
- **figures/O3_productivity/**: `O3_species_energy` ★ **NEW clean fig** (diversity ~ MODIS/PML/tower GPP, monotonic, no DHI) · `L07` tower GPP (16 sites) · `L08` tower GPP vs 4 diversity indices
- **results/O3_productivity/**: `species_energy_stats.csv` (NEW: slopes/quad tests) · `tower_gpp_validation.csv` · `plot_modis_gpp.csv` · `plot_pml_gpp.csv` · `dhi_gpp_triangulation.csv` ⚠ **use MODIS/PML rows only**; DHI row is the excluded artifact
- Script: `scripts/obj23/152_species_energy_clean.py`. Stats: MODIS plot β=+0.52 (p=3.6e-11, r=+0.28); PML β=+0.34 (p=2.6e-5); tower β=+0.26 (n=16, r=+0.21, n.s.). All slopes>0, no inverted-U.

## O4 — Contextual modulation → disturbance, land use, stand age
**Finding**: disturbance loses ~1 rare species (dominance resilient); spectral–alpha coupling strong right after disturbance, weakens with recency; fragmentation → diversity −0.30~−0.35, SAVI×heterogeneity → turnover +0.21; RS–diversity coupling stronger in young stands (age main effect null).
- **figures/O4_context/**: `G03` disturbance map · `G04` coupling in disturbed · `G05` diversity by disturbance · `H01` severity · `H02` recency/recovery · `H03` moderation · `I01` BACI · `I02` severity dose · `I03` change coupling · `K02` land-use heterogeneity · `K03` heterogeneity moderation
- **results/O4_context/**: `obj23_disturbance.csv` · `obj2_heterogeneity.csv` · `obj2_severity_recency.csv` · `baci.csv` · `baci_sensitivity.csv` · `stand_age_moderation_{gami,30m}.csv` · `age_binned_moderation.csv`

---

## EXCLUDED — artifact analyses (deliberately NOT copied)
Per decision (2026-07-21): wrong/artifact analyses are not featured; this is not a caveat paper.
- **DHI greenness "productivity hump"** (greenness artifact, reversed by independent GPP): figures `G01`, `G02`, `L01`, `L02`, `L03`, `L04`, `L05`, `L06`; results `dhi_gpp_validation.csv`, `dhi_hump_confound.csv`, `obj23_dhi_coeff.csv`, `nonlinear.csv`, `gam.csv`.
- **Age–diversity "hump"** (biogeographic artifact, dies under forest-type control): figures `N09`, `N10`; result `foresttype_hump.csv`.

## Provenance
Originals in `NEON_v2/figures/` and `NEON_v2/results/`; superseded per-paper folders `P1_diversity_pattern/`, `P2_structural_process/`, `P3_productivity_hump/` retained for history. Canonical manuscript = `MANUSCRIPT.md` in this folder.
