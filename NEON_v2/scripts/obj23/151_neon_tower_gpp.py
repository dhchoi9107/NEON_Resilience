"""
Validate DHI/MODIS/PML against NEON eddy-covariance TOWER GPP (gold standard).
INPUT: AmeriFlux FLUXNET-1F files for NEON sites, unzipped into  data/ameriflux/
       (annual files: AMF_US-xXX_FLUXNET-1F_*_YY-*.csv, column GPP_NT_VUT_REF).
Aggregates tower GPP per site, compares against the 3 satellite proxies, and
tests site-level diversity ~ tower GPP (hump vs monotonic).
Output: data/site_tower_gpp.csv, results/tower_gpp_validation.csv, figures/L07_tower_gpp.png
"""
import os, glob, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")
AMF=os.path.join(D,"ameriflux")

# AmeriFlux US-x code -> NEON siteID (verify against downloaded filenames; unmatched are reported)
AMF2NEON={'US-xAB':'ABBY','US-xBR':'BART','US-xBL':'BLAN','US-xGR':'GRSM','US-xHA':'HARV',
 'US-xJE':'JERC','US-xML':'MLBS','US-xSB':'OSBS','US-xRM':'RMNP','US-xSC':'SCBI','US-xSE':'SERC',
 'US-xSP':'SOAP','US-xST':'STEI','US-xTA':'TALL','US-xTE':'TEAK','US-xTR':'TREE','US-xUN':'UNDE','US-xWR':'WREF'}

def load_tower_gpp():
    rows=[]
    # annual FLUXNET-1F files end in _YY_ (yearly). GPP in gC m-2 yr-1.
    files=glob.glob(os.path.join(AMF,"**","*FLUXNET*_YY_*.csv"),recursive=True)+\
          glob.glob(os.path.join(AMF,"*FLUXNET*_YY_*.csv"))
    print(f"annual FLUXNET files found: {len(files)}")
    for fp in sorted(set(files)):
        code=next((c for c in AMF2NEON if c in os.path.basename(fp)),None)
        if not code: print("  (unmapped file)",os.path.basename(fp)); continue
        df=pd.read_csv(fp,na_values=[-9999])
        col='GPP_NT_VUT_REF' if 'GPP_NT_VUT_REF' in df else ('GPP_DT_VUT_REF' if 'GPP_DT_VUT_REF' in df else None)
        if not col: print("  (no GPP col)",os.path.basename(fp)); continue
        rows.append((AMF2NEON[code], df[col].mean()))
    return pd.DataFrame(rows,columns=['siteID','tower_gpp']).groupby('siteID',as_index=False).mean()

if not os.path.isdir(AMF) or not glob.glob(os.path.join(AMF,"**","*FLUXNET*"),recursive=True):
    print(f"!! AmeriFlux 데이터 없음. data/ameriflux/ 에 FLUXNET-1F 파일을 풀어주세요 (README 참고). 종료."); raise SystemExit
tg=load_tower_gpp(); tg.to_csv(os.path.join(D,"site_tower_gpp.csv"),index=False)
print(f"사이트 타워 GPP: {len(tg)}개\n{tg.round(1).to_string(index=False)}")

# merge site-level: diversity + 3 satellite proxies + tower GPP
base=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
base=base[base.sample_coverage>=0.9]
for f,c in [("plot_dhi_sentinel.csv","dhi_cum_mean"),("plot_modis_gpp.csv","modis_gpp"),("plot_pml_gpp.csv","pml_gpp")]:
    base=base.merge(pd.read_csv(os.path.join(D,f))[['plotID',c]],on='plotID',how='left')
site=base.groupby('siteID').agg(Hill_q1=('Hill_q1','mean'),DHI=('dhi_cum_mean','median'),
        MODIS=('modis_gpp','median'),PML=('pml_gpp','median')).reset_index().merge(tg,on='siteID')
print(f"\n검증 사이트 수: {len(site)}")

# which satellite proxy best matches tower GPP?
print("\n=== 위성 프록시 vs 타워 GPP 상관 (사이트) ===")
rows=[]
for p in ['DHI','MODIS','PML']:
    r=st.pearsonr(site[p],site.tower_gpp); print(f"  {p:6s}-tower: r={r[0]:+.2f} p={r[1]:.2g}"); rows.append(dict(proxy=p,r_vs_tower=r[0],p=r[1]))
# diversity ~ tower GPP: hump or monotonic?
import statsmodels.formula.api as smf
s=site.copy(); s['zg']=(s.tower_gpp-s.tower_gpp.mean())/s.tower_gpp.std()
m=smf.ols("Hill_q1 ~ zg + I(zg**2)",s).fit()
k=[c for c in m.params.index if '** 2' in c][0]
print(f"\n=== 다양성 ~ 타워 GPP (사이트 n={len(s)}) ===")
print(f"  선형 β={m.params['zg']:+.3f} p={m.pvalues['zg']:.2g} | 2차 β={m.params[k]:+.3f} p={m.pvalues[k]:.2g} -> {'혹형∩' if m.params[k]<0 else '단조/U'}")
rows.append(dict(proxy='tower_quad',r_vs_tower=m.params[k],p=m.pvalues[k]))
pd.DataFrame(rows).to_csv(os.path.join(R,"tower_gpp_validation.csv"),index=False)

# figure
fig,ax=plt.subplots(1,2,figsize=(13,5.2))
ax[0].scatter(site.MODIS,site.tower_gpp,s=50,color='#1565c0',label=f'MODIS r={st.pearsonr(site.MODIS,site.tower_gpp)[0]:.2f}')
ax[0].scatter(site.DHI*site.tower_gpp.mean()/site.DHI.mean(),site.tower_gpp,s=0)  # spacer
for _,r in site.iterrows(): ax[0].annotate(r.siteID,(r.MODIS,r.tower_gpp),fontsize=7)
ax[0].set_xlabel("MODIS GPP"); ax[0].set_ylabel("NEON 타워 GPP"); ax[0].set_title("타워 vs MODIS GPP"); ax[0].legend()
xx=np.linspace(s.tower_gpp.min(),s.tower_gpp.max(),50); zz=(xx-s.tower_gpp.mean())/s.tower_gpp.std()
ax[1].scatter(site.tower_gpp,site.Hill_q1,s=50,color='#c62828')
for _,r in site.iterrows(): ax[1].annotate(r.siteID,(r.tower_gpp,r.Hill_q1),fontsize=7)
ax[1].plot(xx,m.params['Intercept']+m.params['zg']*zz+m.params[k]*zz**2,'k-',lw=2)
ax[1].set_xlabel("NEON 타워 GPP"); ax[1].set_ylabel("Hill q1"); ax[1].set_title("다양성 ~ 타워 GPP (gold standard)")
fig.suptitle("L07. NEON 타워 GPP 검증 — 혹형(∩)인가 단조인가",fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(F,"L07_tower_gpp.png"),dpi=120,bbox_inches='tight'); plt.close()
print("saved L07_tower_gpp.png")
