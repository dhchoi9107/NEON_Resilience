"""
Verify whether the DHI productivity-diversity HUMP is a biogeographic (forest-type)
confound (parallel to the age-hump test in 146).
  - quadratic (hump) p under: raw / +domain / +forest_type / +evergreenness(DHI_min) / +site random
  - within-biome test: does the hump survive inside deciduous (high-div biome) alone?
Output: results/dhi_hump_confound.csv, figures/L04_dhi_hump_confound.png
"""
import os, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")
FT={'ABBY':'evergreen','WREF':'evergreen','TEAK':'evergreen','RMNP':'evergreen','SOAP':'evergreen',
    'JERC':'evergreen','OSBS':'evergreen','TALL':'evergreen','BART':'deciduous','BLAN':'deciduous',
    'SCBI':'deciduous','SERC':'deciduous','MLBS':'deciduous','ORNL':'deciduous','GRSM':'deciduous',
    'HARV':'mixed','UNDE':'mixed','STEI':'mixed','TREE':'mixed'}
ftc={'evergreen':'#1b5e20','deciduous':'#e65100','mixed':'#6a1b9a'}
X='dhi_cum_mean'

df=pd.read_csv(f"{D}/FINAL_v2_pooled.csv").merge(
   pd.read_csv(f"{D}/plot_dhi_sentinel.csv")[['plotID',X,'dhi_min_mean']],on='plotID',how='left')
df=df[df.sample_coverage>=0.9].copy(); df['ft']=df.siteID.map(FT)
df['zx']=(df[X]-df[X].mean())/df[X].std()
DIV={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'turnover','LCBD_nestedness_rare':'nestedness'}

def qp(f,dd,re=True):
    try:
        m=smf.mixedlm(f,dd,groups=dd.siteID).fit(reml=True,method='lbfgs') if re else smf.ols(f,dd).fit()
        for k in m.pvalues.index:
            if 'zx ** 2' in k: return m.pvalues[k],m.params[k]
    except Exception: pass
    return np.nan,np.nan

rows=[]
print("=== DHI 생산성 혹형(2차) p — 통제 강화 순서 ===")
print(f"{'idx':11s} {'raw':>9s} {'+domain':>9s} {'+ftype':>9s} {'+evergrn':>9s} {'+site(RE)':>10s} {'decid내부':>9s}")
for c,l in DIV.items():
    dd=df[[c,'zx','siteID','domain','ft','dhi_min_mean']].dropna()
    pr,br=qp(f"{c} ~ zx + I(zx**2)",dd,False)
    pdm,_=qp(f"{c} ~ zx + I(zx**2) + C(domain)",dd,False)
    pft,_=qp(f"{c} ~ zx + I(zx**2) + C(ft)",dd,False)
    pev,_=qp(f"{c} ~ zx + I(zx**2) + dhi_min_mean",dd,False)
    prs,_=qp(f"{c} ~ zx + I(zx**2) + C(ft)",dd,True)          # forest type + site random
    dec=dd[dd.ft=='deciduous']
    pdc,_=qp(f"{c} ~ zx + I(zx**2)",dec,True) if len(dec)>40 else (np.nan,np.nan)
    rows.append(dict(index=l,raw=pr,domain=pdm,ftype=pft,evergreen=pev,ft_site=prs,decid_only=pdc,n=len(dd),n_decid=len(dec)))
    print(f"{l:11s} {pr:9.2g} {pdm:9.2g} {pft:9.2g} {pev:9.2g} {prs:10.2g} {pdc:9.2g}")
pd.DataFrame(rows).to_csv(f"{R}/dhi_hump_confound.csv",index=False)

# figure L04: DHI vs diversity by forest type, per-type quadratic
fig,axes=plt.subplots(1,2,figsize=(15,6))
for ax,c,l in [(axes[0],'Hill_q1','Hill q1 (알파)'),(axes[1],'LCBD_turnover_rare','LCBD turnover (베타)')]:
    d=df[[c,X,'ft']].dropna()
    for ft in ['evergreen','deciduous','mixed']:
        dd=d[d.ft==ft]
        ax.scatter(dd[X],dd[c],s=14,alpha=.30,color=ftc[ft])
        if len(dd)>10:
            xx=np.linspace(dd[X].min(),dd[X].max(),40); cf=np.polyfit(dd[X],dd[c],2)
            ax.plot(xx,np.polyval(cf,xx),'-',color=ftc[ft],lw=2.6,label=ft)
    xx=np.linspace(d[X].min(),d[X].max(),60); cf=np.polyfit(d[X],d[c],2)
    ax.plot(xx,np.polyval(cf,xx),'k--',lw=2,alpha=.7,label='전체(교란=겉보기 혹형)')
    ax.set_xlabel("DHI cumulative (생산성)"); ax.set_ylabel(c); ax.set_title(l); ax.legend(fontsize=9); ax.grid(alpha=.2)
fig.suptitle("L04. DHI 생산성–다양성 '혹형'은 생물지리 교란 — type 내부엔 혹형 없음(전체 검은점선 vs type별 색)",fontsize=13)
fig.tight_layout(); fig.savefig(f"{F}/L04_dhi_hump_confound.png",dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved L04_dhi_hump_confound.png")
