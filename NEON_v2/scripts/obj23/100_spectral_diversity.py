"""
Obj 1 gap fill — SPECTRAL DIVERSITY (RS-derived diversity), not just VI indices.
Spectral diversity = within-plot spectral heterogeneity: Rao's Q, spectral CV, Shannon,
FRic/FDiv/FEve (PCA functional diversity of pixel spectra). These measure RS-derived diversity;
relate to TAXONOMIC diversity (alpha Hill, beta LCBD turnover/nestedness) at multi scales.
Source: plot_spectral_1m (10 m), filtered to the 8 NEON-BRDF years (2016-18,22-26) for consistency.
Output: FINAL_v2_specdiv.csv + results/obj1_specdiv.csv
"""
import os, numpy as np, pandas as pd, statsmodels.formula.api as smf
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
R=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\results"
def z(s): return (s-s.mean())/s.std()
BRDF_YEARS={2016,2017,2018,2022,2023,2024,2025,2026}
SD=['rao_q','spectral_cv','spectral_shannon','spectral_FRic','spectral_FDiv','spectral_FEve']

sp=pd.read_csv("E:/neon_lidar/spectral_diversity/plot_spectral_1m.csv")
sp=sp[(sp['grain_m']==10)&(sp['year'].isin(BRDF_YEARS))]
div=sp.groupby('plotID')[SD].mean().reset_index()
div.columns=['plotID']+['specdiv_'+c for c in SD]
print(f"spectral diversity (8 BRDF yr): {len(div)} plots | cols {list(div.columns[1:])}")

base=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
m=base.merge(div,on='plotID',how='left')
m.to_csv(os.path.join(D,"FINAL_v2_specdiv.csv"),index=False)
mc=m[m['sample_coverage']>=0.9].copy()
print(f"cov>=0.9 {len(mc)} | with spectral diversity {mc['specdiv_rao_q'].notna().sum()}")

RESP={'Hill_q1':'Hill q1(alpha)','Hill_q2':'Hill q2(alpha)',
      'LCBD_turnover_rare':'LCBD turnover(beta)','LCBD_nestedness_rare':'LCBD nestedness(beta)'}
for r in RESP: mc[f'z_{r}']=z(mc[r])
SDP=['specdiv_'+c for c in SD]
STRUCT_DIV=['FHD_mean','Rumple_mean','Gini_mean','VCI_mean','Vert_SD_mean']  # structural diversity for comparison
rows=[]
for resp,rl in RESP.items():
    for col,kind in [(c,'spectral_div') for c in SDP]+[(c,'structural_div') for c in STRUCT_DIV]:
        d=mc[[f'z_{resp}',col,'siteID','domain']].dropna()
        if len(d)<40 or d['siteID'].nunique()<3: continue
        d=d.copy(); d['zx']=z(d[col])
        try:
            mm=smf.mixedlm(f"z_{resp} ~ zx + C(domain)",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
            if mm.converged:
                rows.append(dict(response=rl,predictor=col.replace('specdiv_','').replace('_mean',''),
                    kind=kind,beta=mm.fe_params['zx'],p=mm.pvalues['zx'],n=len(d)))
        except Exception: pass
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"obj1_specdiv.csv"),index=False)
print("\n=== Obj1: 종다양성 ~ 분광다양성 / 구조다양성 (p<0.05) ===")
for resp,rl in RESP.items():
    d=res[(res.response==rl)&(res.p<0.05)].sort_values('beta',key=abs,ascending=False)
    print(f"\n{rl}:")
    for _,x in d.iterrows():
        print(f"  {x['predictor']:16s} [{x['kind']:13s}] beta={x['beta']:+.3f} p={x['p']:.3g}")
    if d.empty: print("  (유의 없음)")
