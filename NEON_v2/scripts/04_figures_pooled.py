"""
NEON_v2 STEP 4 — Pooled figures (all from FINAL_v2_pooled.csv, single NEON-BRDF VI source).
  F01_forest_<resp>.png   : each response ~ RS (mean|sd|trend) panels  [structural vs VI]
  F02_scatter_<resp>.png  : response vs all RS (mean) scatter grid
  F03_feature_effect.png  : max-R2 by feature (mean/sd/trend) — does interannual SD/trend add?
  F04_variance_decomp.png : domain>site>plot>within(year) variance partition
  F05_by_domain.png / F06_by_site.png : diversity distributions across scales
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=115,bbox_inches='tight'); plt.close(fig); print(" saved",n)
def z(s): return (s-s.mean())/s.std()

df=pd.read_csv(os.path.join(D,"FINAL_v2_pooled.csv"))
df=df[df['sample_coverage']>=0.9].copy()
print("plots(cov>=0.9):",len(df),"sites:",df['siteID'].nunique())
STRUCT=['LAI','FHD','VCI','Canopy_Ht','Max_Ht','Rumple','Rugosity','Deep_Gap','Vert_SD','Vert_CV','Gini','Q95','Ht_Ratio']
VI=['NDVI','EVI','ARVI','SAVI']
PRED=STRUCT+VI; FEAT=['mean','sd','trend']
RESP={'Hill_q1':'Hill q1','Hill_q2':'Hill q2','LCBD_turnover_rare':'LCBD Turnover','LCBD_nestedness_rare':'LCBD Nestedness'}
catc={'Structural':'#c62828','VI':'#2e7d32'}
res=pd.read_csv(os.path.join(R,"v2_coeff.csv"))
rmap={'Hill q1':'Hill q1','Hill q2':'Hill q2','Turnover(rare)':'LCBD Turnover','Nestedness(rare)':'LCBD Nestedness'}
res['response']=res['response'].map(lambda x: rmap.get(x,x))

# ---- F01: forest per response, 3 feature panels, consistent order ----
for resp,rl in RESP.items():
    sub=res[res['response']==rl]
    if sub.empty: continue
    order=sub[sub['feature']=='mean'].sort_values('beta')['predictor'].tolist()
    catmap=sub.drop_duplicates('predictor').set_index('predictor')['category'].to_dict()
    fig,axes=plt.subplots(1,3,figsize=(17,7.5),sharey=True)
    for ax,feat in zip(axes,FEAT):
        d=sub[sub['feature']==feat].set_index('predictor')
        for i,p in enumerate(order):
            if p not in d.index: continue
            b,se,pv=d.loc[p,'beta'],d.loc[p,'se'],d.loc[p,'p']
            ax.errorbar(b,i,xerr=1.96*se,fmt='none',ecolor='gray',zorder=2)
            ax.scatter(b,i,c=catc[catmap[p]],s=70,alpha=1 if pv<0.05 else .25,zorder=3)
        ax.set_yticks(range(len(order)));ax.set_yticklabels(order,fontsize=8)
        ax.axvline(0,color='k',ls='--',lw=1);ax.set_title(f"{rl} ~ RS {feat}");ax.set_xlabel("beta (z)");ax.grid(axis='x',alpha=.3)
    axes[-1].legend(handles=[Patch(color=v,label=k) for k,v in catc.items()],loc='lower right')
    fig.suptitle(f"{rl} ~ per-year RS features (mean / interannual SD / trend)   [VI=NEON .002 BRDF]",fontsize=13)
    fig.tight_layout(); save(fig,f"F01_forest_{resp}.png")

# ---- F02: scatter grids (response vs RS mean) ----
for resp,rl in RESP.items():
    cols_=[f"{c}_mean" for c in PRED]; n=len(cols_); ncol=5; nrow=int(np.ceil(n/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(3.6*ncol,3.0*nrow))
    for ax,c in zip(axes.ravel(),cols_):
        d=df[[c,resp]].dropna(); d=d[np.isfinite(d[c])&np.isfinite(d[resp])]
        if len(d)<10: ax.axis('off'); continue
        col=catc['VI'] if c.replace('_mean','') in VI else catc['Structural']
        ax.scatter(d[c],d[resp],s=8,alpha=0.25,color=col,edgecolors='none')
        sl,ic,r,p,se=st.linregress(d[c],d[resp]); xx=np.linspace(d[c].min(),d[c].max(),40)
        ax.plot(xx,ic+sl*xx,'k-',lw=1.8)
        star='***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
        ax.set_title(f"{c.replace('_mean','')}  r={r:+.2f}{star}",fontsize=9);ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]: ax.axis('off')
    fig.suptitle(f"{rl} vs RS indices (mean)   red=LiDAR structural, green=vegetation index",fontsize=14)
    fig.tight_layout(); save(fig,f"F02_scatter_{resp}.png")

# ---- F03: does interannual SD / trend add beyond mean? max-R2 per feature ----
fig,axes=plt.subplots(1,4,figsize=(20,4.6),sharey=False)
for ax,(resp,rl) in zip(axes,RESP.items()):
    sub=res[res['response']==rl]
    vals=[]; labs=[]
    for feat in FEAT:
        d=sub[sub['feature']==feat]
        if len(d):
            top=d.loc[d['r2'].idxmax()]; vals.append(top['r2']); labs.append(f"{feat}\n{top['predictor']}")
        else: vals.append(0); labs.append(feat)
    bars=ax.bar(range(3),vals,color=['#1565c0','#ef6c00','#6a1b9a'])
    ax.set_xticks(range(3));ax.set_xticklabels(labs,fontsize=9);ax.set_title(rl,fontsize=12)
    ax.set_ylabel("max R²");ax.grid(axis='y',alpha=.3)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v,f"{v:.3f}",ha='center',va='bottom',fontsize=9)
fig.suptitle("F03. Best single-predictor R² by temporal feature — does interannual SD / trend add beyond the multi-year mean?",fontsize=13)
fig.tight_layout(); save(fig,"F03_feature_effect.png")

# ---- F04: hierarchical variance decomposition ----
def var_decomp(data,value,levels):
    d=data[[value]+levels].dropna().copy(); x=d[value].values; N=len(x); grand=x.mean()
    SST=((x-grand)**2).sum(); comp={}; prev=np.full(N,grand)
    for lv in levels:
        gm=d.groupby(levels[:levels.index(lv)+1])[value].transform('mean').values
        comp[lv]=(((gm-prev)**2).sum())/SST if SST>0 else 0; prev=gm
    comp['within']=(((x-prev)**2).sum())/SST if SST>0 else 0
    return comp
rows=[]
for resp,rl in RESP.items():
    c=var_decomp(df,resp,['domain','siteID','plotID']); rows.append(dict(variable=rl,**c))
vd=pd.DataFrame(rows).set_index('variable')
fig,ax=plt.subplots(figsize=(9,5))
cols=['domain','siteID','plotID','within']; colors=['#b71c1c','#ef6c00','#fbc02d','#cfd8dc']
bottom=np.zeros(len(vd))
for c,cl in zip(cols,colors):
    ax.bar(vd.index,vd[c]*100,bottom=bottom*100,label=c,color=cl); bottom+=vd[c].values
ax.set_ylabel("% variance");ax.set_title("F04. Diversity variance partition (domain>site>plot)");ax.legend();ax.set_ylim(0,100)
fig.tight_layout(); save(fig,"F04_variance_decomp.png")
vd.to_csv(os.path.join(R,"v2_variance_decomp.csv"))

# ---- F05/F06: diversity across domain / site ----
for resp,rl in [('Hill_q1','Hill q1'),('Hill_q2','Hill q2')]:
    fig,ax=plt.subplots(figsize=(13,5))
    order=df.groupby('domain')[resp].median().sort_values().index
    data=[df[df['domain']==dm][resp].dropna() for dm in order]
    ax.boxplot(data,labels=order,showfliers=False)
    ax.set_title(f"F05. {rl} by domain");ax.set_ylabel(rl);ax.grid(axis='y',alpha=.3)
    fig.tight_layout(); save(fig,f"F05_{resp}_by_domain.png")
fig,ax=plt.subplots(figsize=(16,5))
order=df.groupby('siteID')['Hill_q1'].median().sort_values().index
ax.boxplot([df[df['siteID']==s]['Hill_q1'].dropna() for s in order],labels=order,showfliers=False)
ax.set_title("F06. Hill q1 by site");ax.set_ylabel("Hill q1");ax.tick_params(axis='x',rotation=45);ax.grid(axis='y',alpha=.3)
fig.tight_layout(); save(fig,"F06_Hill_q1_by_site.png")
print("DONE pooled figures")
