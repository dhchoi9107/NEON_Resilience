"""
Obj 3 — Dynamic Habitat Indices (DHI) from Sentinel-2 intra-annual NDVI (proper seasonal DHI).
Three DHI components (Coops/Radeloff) per plot-year, then plot-level mean + temporal trend.
  cumulative : annual mean greenness (total productivity proxy)
  minimum    : annual minimum (least-productive-period cover)
  variation  : seasonal CV = SD/mean (seasonality)
Robust to irregular sampling: monthly medians, require >=8 valid obs spanning >=5 months.
Reads NEON_v2/data/sentinel/s2_<SITE>.csv -> writes plot_dhi_sentinel.csv (plot-level).
"""
import os, glob, numpy as np, pandas as pd
SDIR=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\sentinel"
OUT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_dhi_sentinel.csv"

files=sorted(glob.glob(os.path.join(SDIR,"s2_*.csv")))
print(f"site files: {len(files)}")
df=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
df=df[(df.ndvi>=-0.1)&(df.ndvi<=1.0)&(df.npix>=1)].copy()
df['month']=pd.to_datetime(df['date']).dt.month
print(f"obs: {len(df)} | plots: {df.plotID.nunique()} | years {df.year.min()}-{df.year.max()}")

# ---- per plot-year DHI ----
py=[]
for (pid,yr),g in df.groupby(['plotID','year']):
    if len(g)<8 or g['month'].nunique()<5: continue
    mm=g.groupby('month')['ndvi'].median()
    cum=float(mm.mean()); mn=float(mm.min()); var=float(mm.std()/mm.mean()) if mm.mean()>0 else np.nan
    py.append((pid,yr,cum,mn,var,len(g)))
pdf=pd.DataFrame(py,columns=['plotID','year','dhi_cum','dhi_min','dhi_var','n_obs'])
print(f"plot-year DHI: {len(pdf)} | plots {pdf.plotID.nunique()}")

# ---- plot-level: mean across years + trend ----
def slope(s,x):
    return np.polyfit(x,s,1)[0] if len(x)>=3 else np.nan
rows=[]
for pid,g in pdf.groupby('plotID'):
    g=g.sort_values('year'); yr=g['year'].values
    rec={'plotID':pid,'dhi_nyears':len(g)}
    for c in ['dhi_cum','dhi_min','dhi_var']:
        rec[f'{c}_mean']=float(g[c].mean())
        rec[f'{c}_trend']=float(slope(g[c].values,yr))
    rows.append(rec)
out=pd.DataFrame(rows)
out.to_csv(OUT,index=False)
print(f"DONE -> {OUT}\nplot-level DHI: {out.shape}, plots {out.plotID.nunique()}")
print(out[['dhi_cum_mean','dhi_min_mean','dhi_var_mean']].describe().round(3).to_string())
