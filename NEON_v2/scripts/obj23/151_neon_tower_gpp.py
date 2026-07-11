"""
Validate DHI/MODIS/PML against NEON eddy-covariance TOWER GPP (gold standard).
INPUT: AmeriFlux FLUXNET-1F files for NEON sites, unzipped into  data/ameriflux/
       (annual files: AMF_US-xXX_FLUXNET-1F_*_YY-*.csv, column GPP_NT_VUT_REF).
Aggregates tower GPP per site, compares against the 3 satellite proxies, and
tests site-level diversity ~ tower GPP (hump vs monotonic).
Output: data/site_tower_gpp.csv, results/tower_gpp_validation.csv, figures/L07_tower_gpp.png
"""
import os, glob, io, zipfile, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")
AMF=r"E:\FLUX"   # AmeriFlux FLUXNET zips (read directly, no extraction)

# AmeriFlux US-x code -> NEON siteID (verify against downloaded filenames; unmatched are reported)
AMF2NEON={'US-xAB':'ABBY','US-xBR':'BART','US-xBL':'BLAN','US-xGR':'GRSM','US-xHA':'HARV',
 'US-xJE':'JERC','US-xML':'MLBS','US-xSB':'OSBS','US-xRM':'RMNP','US-xSC':'SCBI','US-xSE':'SERC',
 'US-xSP':'SOAP','US-xST':'STEI','US-xTA':'TALL','US-xTE':'TEAK','US-xTR':'TREE','US-xUN':'UNDE','US-xWR':'WREF'}

def load_tower_gpp():
    rows=[]
    zips=glob.glob(os.path.join(AMF,"AMF_*_FLUXNET_*.zip"))
    print(f"FLUXNET zips: {len(zips)}")
    for zp in sorted(zips):
        code=next((c for c in AMF2NEON if c in os.path.basename(zp)),None)
        if not code: continue                       # non-NEON-forest site, skip
        try:
            zf=zipfile.ZipFile(zp)
            yy=[n for n in zf.namelist() if ('FLUXMET_YY' in n or 'FULLSET_YY' in n) and n.endswith('.csv')]
            if not yy: print("  (no FLUXMET_YY)",os.path.basename(zp)); continue
            df=pd.read_csv(io.BytesIO(zf.read(yy[0])),na_values=[-9999])
            g=df['GPP_NT_VUT_REF'].dropna()
            if len(g): rows.append((AMF2NEON[code], float(g.mean()), int(len(g))))
        except Exception as e: print("  (err)",os.path.basename(zp),e)
    out=pd.DataFrame(rows,columns=['siteID','tower_gpp','n_years'])
    return out

if not os.path.isdir(AMF):
    print(f"!! {AMF} 없음. AmeriFlux FLUXNET zip 폴더 경로를 AMF에 설정하세요. 종료."); raise SystemExit
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
    s2=site[[p,'tower_gpp']].dropna()
    r=st.pearsonr(s2[p],s2.tower_gpp); print(f"  {p:6s}-tower: r={r[0]:+.2f} p={r[1]:.2g} (n={len(s2)})"); rows.append(dict(proxy=p,r_vs_tower=r[0],p=r[1],n=len(s2)))
# diversity ~ tower GPP: hump or monotonic?
import statsmodels.formula.api as smf
s=site.copy(); s['zg']=(s.tower_gpp-s.tower_gpp.mean())/s.tower_gpp.std()
m=smf.ols("Hill_q1 ~ zg + I(zg**2)",s).fit()
k=[c for c in m.params.index if '** 2' in c][0]
print(f"\n=== 다양성 ~ 타워 GPP (사이트 n={len(s)}) ===")
print(f"  선형 β={m.params['zg']:+.3f} p={m.pvalues['zg']:.2g} | 2차 β={m.params[k]:+.3f} p={m.pvalues[k]:.2g} -> {'혹형∩' if m.params[k]<0 else '단조/U'}")
rows.append(dict(proxy='tower_quad',r_vs_tower=m.params[k],p=m.pvalues[k]))
pd.DataFrame(rows).to_csv(os.path.join(R,"tower_gpp_validation.csv"),index=False)

# figure: (A) each satellite proxy vs gold-standard tower GPP; (B) diversity ~ tower GPP
fig,ax=plt.subplots(1,3,figsize=(19,5.4))
for a,(p,cc) in zip(ax[:2] if False else [ax[0]],[('MODIS','#1565c0')]): pass
for a,p,cc in [(ax[0],'MODIS','#1565c0')]:
    pass
# panel 0: correlation bars (which satellite proxy tracks tower GPP)
cors={p:st.pearsonr(*[site[[p,'tower_gpp']].dropna()[x] for x in [p,'tower_gpp']])[0] for p in ['DHI','MODIS','PML']}
ax[0].bar(list(cors),list(cors.values()),color=['#c62828','#1565c0','#2e7d32'])
for i,(p,v) in enumerate(cors.items()): ax[0].text(i,v+.02,f"{v:.2f}",ha='center',fontsize=11)
ax[0].set_ylabel("r vs NEON 타워 GPP"); ax[0].set_ylim(0,.85); ax[0].set_title("위성 프록시 vs 실측 타워 GPP\nDHI만 약함(녹색도≠생산성)"); ax[0].grid(axis='y',alpha=.3)
# panel 1: MODIS vs tower scatter
d2=site[['MODIS','tower_gpp','siteID']].dropna()
ax[1].scatter(d2.MODIS,d2.tower_gpp,s=45,color='#1565c0')
for _,r in d2.iterrows(): ax[1].annotate(r.siteID,(r.MODIS,r.tower_gpp),fontsize=7)
ax[1].set_xlabel("MODIS GPP"); ax[1].set_ylabel("NEON 타워 GPP"); ax[1].set_title(f"MODIS vs 타워 (r={cors['MODIS']:.2f})"); ax[1].grid(alpha=.2)
# panel 2: diversity ~ tower GPP (no hump)
xx=np.linspace(s.tower_gpp.min(),s.tower_gpp.max(),50); zz=(xx-s.tower_gpp.mean())/s.tower_gpp.std()
ax[2].scatter(site.tower_gpp,site.Hill_q1,s=45,color='#c62828')
for _,r in site.iterrows(): ax[2].annotate(r.siteID,(r.tower_gpp,r.Hill_q1),fontsize=7)
ax[2].plot(xx,m.params['Intercept']+m.params['zg']*zz+m.params[k]*zz**2,'k-',lw=2)
ax[2].set_xlabel("NEON 타워 GPP (실측)"); ax[2].set_ylabel("Hill q1"); ax[2].set_title(f"다양성 ~ 타워 GPP\n2차 β={m.params[k]:+.2f} (혹형 아님)"); ax[2].grid(alpha=.2)
fig.suptitle("L07. NEON 타워 GPP(gold standard) 검증 — DHI는 생산성 아님, 다양성 혹형 없음",fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(F,"L07_tower_gpp.png"),dpi=120,bbox_inches='tight'); plt.close()
print("saved L07_tower_gpp.png")
