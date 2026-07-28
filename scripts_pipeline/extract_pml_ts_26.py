# PML-V2 GPP 2000-2020 annual time series for 26 sites -> mean, trend, interannual SD.
import ee, numpy as np, pandas as pd, os
ee.Initialize(project='geedankook')
D="NEON_v2/data"
pl=pd.read_csv(f"{D}/plot_lonlat_26.csv")
plots=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon,r.lat]),{'plotID':r.plotID}) for r in pl.itertuples()])
PML=ee.ImageCollection("projects/pml_evapotranspiration/PML/OUTPUT/PML_V22a").select('GPP')  # 8-day gC/m2/d, 500m
rows=[]
for y in range(2000,2025):
    img=PML.filterDate(f"{y}-01-01",f"{y}-12-31").mean()  # annual mean 8-day GPP
    fc=img.reduceRegions(plots,ee.Reducer.mean(),500).getInfo()
    for f in fc['features']:
        p=f['properties']; v=p.get('mean',p.get('GPP'))
        if v is not None: rows.append((p['plotID'],y,v))
    print(f"  {y}: done", flush=True)
ts=pd.DataFrame(rows,columns=['plotID','year','pml']).dropna()
ts.to_csv(f"{D}/plot_pml_gpp_annual_26.csv",index=False)
print(f"per-year rows: {len(ts)} | plots: {ts.plotID.nunique()} | years {ts.year.min()}-{ts.year.max()}", flush=True)
# per-plot mean, trend (slope/yr), interannual SD
def agg(g):
    g=g.sort_values('year'); yrs=g.year.values.astype(float); v=g.pml.values.astype(float)
    slope=np.polyfit(yrs,v,1)[0] if len(v)>=5 else np.nan
    return pd.Series({'pml_gpp':v.mean(),'pml_gpp_trend':slope,'pml_gpp_sd':v.std()})
out=ts.groupby('plotID').apply(agg).reset_index()
out.to_csv(f"{D}/plot_pml_gpp_ts_26.csv",index=False)
print(f"\nplot_pml_gpp_ts_26.csv: {len(out)} plots", flush=True)
print("  pml_gpp(2000-2020 mean): 범위 %.2f~%.2f"%(out.pml_gpp.min(),out.pml_gpp.max()), flush=True)
print("  pml_gpp_trend (slope/yr): 중앙 %.4f, 음수 %d / 양수 %d"%(out.pml_gpp_trend.median(),(out.pml_gpp_trend<0).sum(),(out.pml_gpp_trend>0).sum()), flush=True)
