"""
Obj 2 + Obj 3 integrated analysis.
Inputs: FINAL_v2_pooled.csv (diversity + LiDAR + VI), plot_dhi_sentinel.csv (Obj3),
        plot_hansen_loss.csv + plot_mtbs_fire.csv + plot_disturb_s2.csv (Obj2 disturbance).

Obj 3: pooled diversity ~ Sentinel DHI components (cumulative/minimum/variation; mean + trend).
Obj 2: unified disturbance (fire/loss/spectral) -> does disturbance MODERATE the
       RS<->taxonomic relationship? Interaction  diversity ~ z(RS) * disturbed.
Outputs: results/obj23_dhi_coeff.csv, results/obj23_disturbance.csv, merged dataset.
"""
import os, numpy as np, pandas as pd
import statsmodels.formula.api as smf
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"results")
def z(s): return (s-s.mean())/s.std()

df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
dhi=pd.read_csv(os.path.join(D,"plot_dhi_sentinel.csv"))
han=pd.read_csv(os.path.join(D,"plot_hansen_loss.csv"))
mtbs=pd.read_csv(os.path.join(D,"plot_mtbs_fire.csv"))
rob=pd.read_csv(os.path.join(D,"plot_disturbance_robust.csv"))  # ROBUST detection (sev+recency)
m=df.merge(dhi,on='plotID',how='left').merge(han,on='plotID',how='left')\
     .merge(mtbs,on='plotID',how='left').merge(rob[['plotID','sev_spec','sev_loss','severity',
        'd_spec','d_loss','d_fire','disturbed','dist_year','recency']],on='plotID',how='left')

# ---- unified disturbance (robust) ----
m['dist_fire']=m['d_fire'].fillna(0).astype(int)
m['dist_loss']=m['d_loss'].fillna(0).astype(int)
m['dist_spec']=m['d_spec'].fillna(0).astype(int)
m['disturbed']=m['disturbed'].fillna(0).astype(int)
m['dist_type']=np.where(m['dist_fire']==1,'fire',np.where(m['dist_loss']==1,'loss',
               np.where(m['dist_spec']==1,'spectral','none')))
m['years_since_dist']=m['recency']
m.to_csv(os.path.join(D,"FINAL_v2_obj23.csv"),index=False)
mc=m[m['sample_coverage']>=0.9].copy()
print(f"plots(cov>=0.9): {len(mc)} | disturbed: {mc['disturbed'].sum()} "
      f"(fire {mc['dist_fire'].sum()}, loss {mc['dist_loss'].sum()}, spec {mc['dist_spec'].sum()}) "
      f"| with DHI: {mc['dhi_cum_mean'].notna().sum()}")

RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'Turnover','LCBD_nestedness_rare':'Nestedness'}
for r in RESP: mc[f"z_{r}"]=z(mc[r])

def run(resp,formula,d):
    d=d.dropna(subset=[resp]+[c for c in d.columns if c in formula and c!='C(domain)'])
    if len(d)<40 or d['siteID'].nunique()<3: return None
    try: mm=smf.mixedlm(f"{resp} ~ {formula}",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
    except Exception: return None
    return mm

# ===== Obj 3: diversity ~ DHI =====
DHI=['dhi_cum_mean','dhi_min_mean','dhi_var_mean','dhi_cum_trend','dhi_var_trend']
rows=[]
for resp,rl in RESP.items():
    for c in DHI:
        d=mc[[f"z_{resp}",c,'siteID','domain']].dropna()
        if len(d)<40: continue
        d=d.copy(); d['zx']=z(d[c])
        mm=run(f"z_{resp}","zx + C(domain)",d)
        if mm is not None and mm.converged:
            rows.append(dict(response=rl,dhi=c,beta=mm.fe_params['zx'],
                             se=mm.bse_fe['zx'],p=mm.pvalues['zx'],n=len(d)))
obj3=pd.DataFrame(rows); obj3.to_csv(os.path.join(R,"obj23_dhi_coeff.csv"),index=False)
print("\n=== Obj3: diversity ~ Sentinel DHI (p<0.05) ===")
for _,x in obj3[obj3.p<0.05].sort_values('p').iterrows():
    print(f"  {x['response']:10s} ~ {x['dhi']:16s} beta={x['beta']:+.3f} p={x['p']:.3g} n={x['n']}")

# ===== Obj 2: does disturbance moderate RS<->taxonomic coupling? =====
# interaction: diversity ~ z(RS) * disturbed
KEYRS=['SAVI_mean','EVI_mean','Deep_Gap_mean','FHD_mean','LAI_mean','dhi_cum_mean']
rows=[]
for resp,rl in RESP.items():
    for rs in KEYRS:
        d=mc[[f"z_{resp}",rs,'disturbed','siteID','domain']].dropna()
        if len(d)<60 or d['disturbed'].sum()<10: continue
        d=d.copy(); d['zx']=z(d[rs])
        mm=run(f"z_{resp}","zx * disturbed + C(domain)",d)
        if mm is None or not mm.converged: continue
        ik=[k for k in mm.fe_params.index if 'zx:disturbed' in k or 'zx:' in k and 'disturbed' in k]
        inter=mm.fe_params.get('zx:disturbed',np.nan); pint=mm.pvalues.get('zx:disturbed',np.nan)
        rows.append(dict(response=rl,predictor=rs,slope_undist=mm.fe_params['zx'],
                         interaction=inter,p_interaction=pint,n=len(d),n_dist=int(d['disturbed'].sum())))
obj2=pd.DataFrame(rows); obj2.to_csv(os.path.join(R,"obj23_disturbance.csv"),index=False)
print("\n=== Obj2: disturbance x RS interaction on diversity ===")
print("  (interaction p<0.1 = disturbance changes the RS-taxonomic slope)")
for _,x in obj2.sort_values('p_interaction').iterrows():
    sig='*' if x['p_interaction']<0.05 else '.' if x['p_interaction']<0.1 else ' '
    print(f"  {x['response']:10s} ~ {x['predictor']:14s} slope={x['slope_undist']:+.2f} "
          f"interaction={x['interaction']:+.2f}{sig} (p={x['p_interaction']:.2g}) n_dist={x['n_dist']}")
print("\nDONE")
