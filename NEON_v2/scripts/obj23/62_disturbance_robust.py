"""
Obj 2 — ROBUST disturbance severity + recency (replaces noisy -0.08 binary).
Spectral disturbance (from annual NBR series) flagged only if the largest NBR drop is:
  (a) absolutely large (> 0.10), AND
  (b) anomalous vs the plot's OWN year-to-year variability (> 3*MAD), AND
  (c) SUSTAINED (post-drop mean NBR < pre-drop mean - 0.04; not a 1-yr noise dip).
Severity = drop magnitude (continuous). Combined with Hansen loss fraction + MTBS fire.
Recency = years since most recent disturbance (from any source).
Output: plot_disturbance_robust.csv
"""
import numpy as np, pandas as pd
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
ts=pd.read_csv(D+r"\plot_nbr_annual.csv")
han=pd.read_csv(D+r"\plot_hansen_loss.csv")
mtbs=pd.read_csv(D+r"\plot_mtbs_fire.csv")

def detect(g):
    g=g.sort_values('year'); y=g['year'].values; v=g['nbr'].values
    if len(g)<5: return pd.Series(dict(sev_spec=0.0,spec_year=np.nan))
    diffs=np.diff(v); mad=np.median(np.abs(diffs-np.median(diffs)))*1.4826
    i=int(np.argmin(diffs)); drop=diffs[i]; ydrop=int(y[i+1])
    thresh=max(0.10,3*mad)
    pre=v[y<ydrop].mean(); post=v[y>=ydrop].mean()
    sustained=post < pre-0.04
    flag=(drop < -thresh) and sustained
    return pd.Series(dict(sev_spec=float(-drop) if flag else 0.0,
                          spec_year=ydrop if flag else np.nan))
spec=ts.groupby('plotID').apply(detect).reset_index()
print("spectral robust disturbance flagged:",(spec.sev_spec>0).sum(),"/",len(spec),
      "| sev_spec>0 중앙:",round(spec[spec.sev_spec>0].sev_spec.median(),3))

m=spec.merge(han[['plotID','hansen_loss','hansen_lossyear','hansen_loss_frac']],on='plotID',how='outer')\
      .merge(mtbs[['plotID','fire','fire_year']],on='plotID',how='outer')
m['sev_loss']=m['hansen_loss_frac'].fillna(0.0)
m['sev_spec']=m['sev_spec'].fillna(0.0)
# unified: disturbed if any robust source
m['d_spec']=(m['sev_spec']>0).astype(int)
m['d_loss']=(m['hansen_loss'].fillna(0)>0).astype(int)
m['d_fire']=(m['fire'].fillna(0)>0).astype(int)
m['disturbed']=((m.d_spec+m.d_loss+m.d_fire)>0).astype(int)
# most-recent disturbance year across sources
yrs=m[['spec_year','hansen_lossyear','fire_year']].replace(0,np.nan)
m['dist_year']=yrs.max(axis=1)
m['recency']=2025-m['dist_year']               # years since (NaN if undisturbed)
# combined severity z-blend (spectral drop + loss fraction, each scaled 0-1-ish)
m['severity']=m['sev_spec']+m['sev_loss']       # both ~0-1 scale, additive intensity
m.to_csv(D+r"\plot_disturbance_robust.csv",index=False)
print(f"\nDONE plots {len(m)} | disturbed {m.disturbed.sum()} "
      f"(fire {m.d_fire.sum()}, loss {m.d_loss.sum()}, spec {m.d_spec.sum()})")
print("recency(yr) 분포:",m['recency'].describe()[['min','25%','50%','75%','max']].round(1).to_dict())
print("severity 분포(교란만):",m[m.disturbed==1]['severity'].describe()[['min','50%','max']].round(3).to_dict())
