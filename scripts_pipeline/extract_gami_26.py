# GAMI v3.1 forest age (2020, ensemble-mean) for all 26-site plots. Remote zarr (GFZ S3).
# Validated against original plot_stand_age_gami.csv (r=0.83). Re-extracts ALL 26 for consistency.
import xarray as xr, numpy as np, pandas as pd, os
URL="https://s3.gfz-potsdam.de/dog.atlaseo-glm.eo-gridded-data/collections/GAMI/GAMI_v3.1.zarr"
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
pl=pd.read_csv(os.path.join(D,"plot_lonlat_26.csv"))
print(f"plots: {len(pl)} | sites: {pl.siteID.nunique()}", flush=True)
ds=xr.open_zarr(URL, decode_times=True, consolidated=True)
fa=ds['forest_age'].isel(time=1)  # 2020
lats=xr.DataArray(pl.lat.values,dims='p'); lons=xr.DataArray(pl.lon.values,dims='p')
sub=fa.sel(latitude=lats,longitude=lons,method='nearest')   # (members, p)
print("sampling ensemble-mean over 20 members ...", flush=True)
age=sub.mean('members').compute().values
pl['stand_age_gami']=np.round(age,1)
out=os.path.join(D,"plot_stand_age_gami_26.csv")
pl[['plotID','siteID','lon','lat','stand_age_gami']].to_csv(out,index=False)
new=['DELA','LENO','UKFS','YELL','BONA','DEJU','HEAL']
print(f"DONE -> {out}", flush=True)
print(f"non-null: {pl.stand_age_gami.notna().sum()}/{len(pl)} | new-site cov: {pl[pl.siteID.isin(new)].dropna(subset=['stand_age_gami']).siteID.nunique()}/7", flush=True)
print(pl.groupby('siteID').stand_age_gami.median().round(0).to_string())
