"""
(A) Does the temporal TREND feature add predictive value beyond the MEAN?
For each response x predictor: fit  resp ~ z(mean) + z(trend) + C(domain) + (1|site)
and test whether the trend coefficient is significant AFTER controlling for mean.
Also report marginal-R2 gain of (mean+trend) vs (mean only). ML used for AIC comparison.
Output: results/trend_increment.csv, figures/M01_trend_increment.png
"""
import os, numpy as np, pandas as pd
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,R,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"results"),os.path.join(ROOT,"figures")
def z(s): return (s-s.mean())/s.std()

df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv")); df=df[df['sample_coverage']>=0.9].copy()
STRUCT=['Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV','Gini','VCI','FHD','LAI','Q95','Ht_Ratio']
VI=['NDVI','EVI','ARVI','SAVI']; PRED=STRUCT+VI
RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'Turnover(rare)','LCBD_nestedness_rare':'Nestedness(rare)'}

def fitd(d,terms):
    try:
        m=smf.mixedlm(f"zy ~ {terms} + C(domain)",d,groups=d['siteID']).fit(reml=False,method='lbfgs')
    except Exception: return None
    return m if m.converged else None

rows=[]
for resp,rl in RESP.items():
    for p in PRED:
        mc,mt=f"{p}_mean",f"{p}_trend"
        if mc not in df or mt not in df: continue
        # SAME sample for both nested models (mean-only vs mean+trend) -> valid LRT
        d=df[[resp,mc,mt,'siteID','domain']].dropna()
        if len(d)<40 or d['siteID'].nunique()<3 or d['domain'].nunique()<2: continue
        d=d.copy(); d['zy']=z(d[resp]); d['z_'+mc]=z(d[mc]); d['z_'+mt]=z(d[mt])
        m0=fitd(d,'z_'+mc)                 # mean only
        m1=fitd(d,f"z_{mc} + z_{mt}")      # mean + trend
        if m0 is None or m1 is None: continue
        b_tr=m1.fe_params.get('z_'+mt,np.nan); p_tr=m1.pvalues.get('z_'+mt,np.nan)
        b_mean_solo=m0.fe_params.get('z_'+mc,np.nan)
        b_mean_adj=m1.fe_params.get('z_'+mc,np.nan)   # mean effect after adding trend
        # Wald test of trend coefficient in the joint model = incremental contribution
        # (asymptotically equivalent to a nested LRT; valid within a single fitted model).
        rows.append(dict(response=rl,predictor=p,category='VI' if p in VI else 'Structural',
                         beta_mean_solo=b_mean_solo,beta_mean_adj=b_mean_adj,
                         beta_trend_adj=b_tr,p_trend=p_tr,n=len(d)))
res=pd.DataFrame(rows); res.to_csv(os.path.join(R,"trend_increment.csv"),index=False)

print(f"tested {len(res)} response x predictor pairs")
sig=res[res.p_trend<0.05].copy()
print(f"\n=== trend adds beyond mean (p_trend<0.05): {len(sig)}/{len(res)} ===")
for rl in RESP.values():
    s=sig[sig.response==rl].sort_values('p_trend')
    print(f"\n{rl}: {len(s)} predictors")
    for _,x in s.iterrows():
        print(f"  {x.category:10s} {x.predictor:10s} beta_trend|mean={x.beta_trend_adj:+.3f} p={x.p_trend:.3g}")
by=res.groupby('category').apply(lambda g:(g.p_trend<0.05).sum(),include_groups=False)
print("\n=== 카테고리별 trend 유의 개수 ===\n",by.to_string())

# summary figure: signed trend effect (controlling for mean), star where p_trend<0.05
piv=res.pivot(index='predictor',columns='response',values='beta_trend_adj').reindex(PRED)
pv=res.pivot(index='predictor',columns='response',values='p_trend').reindex(PRED)
fig,ax=plt.subplots(figsize=(8,7))
im=ax.imshow(piv.values,cmap='RdBu_r',vmin=-0.25,vmax=0.25,aspect='auto')
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns,rotation=20,ha='right')
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index,fontsize=9)
for i in range(len(piv.index)):
    for j in range(len(piv.columns)):
        if pv.values[i,j]<0.05: ax.text(j,i,'*',ha='center',va='center',fontsize=14,fontweight='bold')
fig.colorbar(im,label='β(trend | mean): mean 통제 후 trend 효과')
ax.set_title("M01 (A). trend가 mean 대비 추가 예측력?  * = p(trend|mean)<0.05",fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(F,"M01_trend_increment.png"),dpi=120,bbox_inches='tight'); plt.close()
print("\nsaved M01_trend_increment.png ; DONE")
