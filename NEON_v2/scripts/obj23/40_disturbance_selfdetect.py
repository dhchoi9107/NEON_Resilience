"""
Obj 2 (self-detect) — disturbance from Sentinel-2 NBR time series (LandTrendr-lite).
Per plot: annual growing-season NBR -> largest year-over-year drop = candidate disturbance.
Captures events any cause (fire/harvest/insect/windthrow) that depress canopy NBR.
Output: plot_disturb_s2.csv (plotID, s2_dist, s2_dist_year, s2_dist_mag, s2_recovery, nbr_trend)
"""
import os, glob, numpy as np, pandas as pd
SDIR=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\sentinel"
OUT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_disturb_s2.csv"
DROP=-0.08   # NBR drop threshold flagging disturbance

df=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(SDIR,"s2_*.csv")))],ignore_index=True)
df=df[(df.nbr>=-1)&(df.nbr<=1)].copy()
df['month']=pd.to_datetime(df['date']).dt.month
gs=df[(df.month>=6)&(df.month<=9)]   # growing-season (peak canopy)
print(f"obs {len(df)} | plots {df.plotID.nunique()}")

rows=[]
for pid,g in gs.groupby('plotID'):
    ay=g.groupby('year')['nbr'].median().sort_index()
    if len(ay)<4:
        rows.append((pid,0,0,np.nan,np.nan,np.nan)); continue
    yrs=ay.index.values; vals=ay.values
    d=np.diff(vals)                     # year-over-year change
    i=int(np.argmin(d)); mag=float(d[i]); dyear=int(yrs[i+1])
    # recovery: slope of NBR after the drop year
    post=ay[ay.index>=dyear]
    rec=float(np.polyfit(post.index,post.values,1)[0]) if len(post)>=3 else np.nan
    trend=float(np.polyfit(yrs,vals,1)[0])
    rows.append((pid,int(mag<DROP),dyear if mag<DROP else 0,round(mag,3),
                 round(rec,4) if rec==rec else np.nan,round(trend,4)))
out=pd.DataFrame(rows,columns=['plotID','s2_dist','s2_dist_year','s2_dist_mag','s2_recovery','nbr_trend'])
out.to_csv(OUT,index=False)
print(f"DONE -> {OUT}\nplots {len(out)} | disturbance flagged: {out.s2_dist.sum()}")
print("flagged year dist:"); print(out[out.s2_dist==1]['s2_dist_year'].value_counts().sort_index().to_string())
