# RS variation hypothesis for beta diversity (within-site)

Question: do plots more dissimilar in remote-sensing (structural/spectral) space also
differ more in species composition? Tested two ways, within-site, biogeography controlled.

## Within-site restricted-permutation Mantel (999 perms, 23 sites, 8,813 plot pairs)
Note: 9% of within-site pairs share no species (Bray-Curtis=1, ceiling) -> rank-based (Spearman) Mantel used.
| RS dissimilarity | vs compositional (Bray–Curtis) | p |
|---|---|---|
| Structural distance (Euclidean, LiDAR metrics) | Spearman Mantel r = +0.32 (Pearson +0.26) | 0.001 |
| Spectral distance (Euclidean, VIs) | Spearman Mantel r = +0.30 (Pearson +0.29) | 0.001 |

## LCBD (plot uniqueness) correlations
| RS LCBD | vs turnover LCBD | vs nestedness LCBD |
|---|---|---|
| structural | r = +0.22 (p<0.001) | +0.02 (n.s.) |
| spectral | r = +0.24 (p<0.001) | −0.10 |

Conclusion: remote-sensing heterogeneity DOES track compositional turnover (modest,
~7% variance, robust p=0.001), but NOT richness difference (nestedness). The signal
appears only when RS is analysed as dissimilarity, not as plot-level predictor levels —
the plot-level framework (§3.1) understates the beta signal because it asks a different
question. Scripts: rs_beta_mantel.py, fig_rs_variation.py -> figures/F3_rs_variation_beta.png

## Saturation model (GDM-family): logit-link GLM, structure + spectral jointly
Handles the Bray-Curtis [0,1] ceiling via logit link; structure & spectral controlled together;
within-site restricted permutation (999) for p.
| predictor (joint) | coefficient | within-site perm p |
|---|---|---|
| structural distance | +0.23 | 0.001 |
| spectral distance   | +0.30 | 0.001 |
| pseudo-R² | 0.114 | |
Both remain independently significant with the ceiling modeled and mutual control — structure
and spectral heterogeneity make non-redundant contributions to compositional turnover.
Script: rs_beta_gdm.py
