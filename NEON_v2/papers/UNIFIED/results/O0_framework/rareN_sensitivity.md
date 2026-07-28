# Supplementary — Rarefaction sensitivity (beta diversity)

Rarefied LCBD recomputed for all 26 sites; within-site Baselga; Monte-Carlo 999 draws (seed 11).
Nested cluster-robust test (domain FE + site-clustered SE), dynamics & productivity block increments.

| RARE_N | n | dynamics→nestedness (dR², p) | productivity→turnover (dR², p) |
|-------:|--:|:--|:--|
| 10 | 579 | 0.056, p=0.006 | 0.015, p=0.006 |
| 15 | 522 | 0.054, p=0.007 | 0.022, p=0.002 |
| 20 | 464 | 0.051, p=0.115 | 0.025, p=0.021 |

Conclusions: **productivity→turnover increment is robust across rarefaction depth** (significant at
RARE_N 10/15/20). The **dynamics→nestedness effect size is stable** (largest unique block at all depths,
dR²=0.051–0.056) but its **statistical significance is power-dependent** — it drops out at RARE_N=20 as
plot retention falls (579→464). Standard adopted: RARE_N=10, N_DRAWS=999.
Monte-Carlo 30→999 stabilization moved the marginal nestedness–productivity wild-bootstrap p from 0.098
to 0.183 (now clearly n.s.); all other conclusions unchanged.
