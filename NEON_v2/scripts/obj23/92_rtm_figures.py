"""
Updated RTM figures: (J05) which structural metrics survive RTM control (impact vs control
slopes) + interaction summary; (J06) worked calculation walkthrough on a real plot.
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures")
def z(s): return (s-s.mean())/s.std()
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=120,bbox_inches='tight'); plt.close(fig); print(" saved",n)

neon=pd.read_csv(os.path.join(D,"plot_disturbance_neon.csv"))
imp=neon[neon.disturbed==1].dropna(subset=['dist_year'])
impmap=imp.set_index('plotID')['dist_year'].to_dict(); sev=neon.set_index('plotID')['severity'].to_dict()
imp['siteID']=imp.plotID.str[:4]; smed=imp.groupby('siteID')['dist_year'].median().to_dict(); gmed=float(imp['dist_year'].median())
def split(pid): return impmap.get(pid,smed.get(pid[:4],gmed))
py=pd.read_csv(os.path.join(D,"per_year_v2.csv"))
METS=['Canopy_Ht','Max_Ht','LAI','Deep_Gap','Gini','FHD','VCI','Vert_SD']
rows=[]
for pid,g in py.groupby('plotID'):
    sy=split(pid); b=g[g.year<sy]; a=g[g.year>=sy]
    if len(b)<1 or len(a)<1: continue
    rec={'plotID':pid,'impact':int(pid in impmap),'severity':sev.get(pid,0) or 0,'pre_complex':b['Rugosity'].mean()}
    for m in METS: rec['d_'+m]=a[m].mean()-b[m].mean()
    rows.append(rec)
ch=pd.DataFrame(rows); ch['siteID']=ch.plotID.str[:4]; ch['zpre']=z(ch['pre_complex'])

def inter_p(m):
    d=ch.dropna(subset=[m])
    try:
        mm=smf.mixedlm(f"{m} ~ impact*zpre + severity",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
        return mm.fe_params['impact:zpre'],mm.pvalues['impact:zpre']
    except: return np.nan,np.nan

# ===== J05: survivors vs RTM =====
show=[('d_Canopy_Ht','ΔCanopy height','RTM 아티팩트'),
      ('d_Deep_Gap','ΔDeep Gap (갭)','RTM 통과 (O)'),
      ('d_Gini','ΔGini (이질성)','RTM 통과 (O)')]
fig,axes=plt.subplots(1,4,figsize=(21,5))
for ax,(m,lab,tag) in zip(axes,show):
    for gflag,c,nm in [(0,'#1565c0','control (RTM baseline)'),(1,'#d32f2f','impact (disturbed)')]:
        d=ch[ch.impact==gflag].dropna(subset=[m,'pre_complex'])
        ax.scatter(d['pre_complex'],d[m],s=14,alpha=.30,color=c)
        sl,ic,r,p,se=st.linregress(d['pre_complex'],d[m]); xx=np.linspace(d.pre_complex.min(),d.pre_complex.max(),20)
        ax.plot(xx,ic+sl*xx,'-',color=c,lw=2.6,label=f"{nm}: slope {sl:+.3f}")
    bi,pi=inter_p(m)
    ax.axhline(0,color='gray',ls=':'); ax.set_xlabel("초기 복잡도 (pre rugosity)"); ax.set_ylabel(lab)
    ax.legend(fontsize=8,loc='best'); ax.set_title(f"{lab}\n{tag}  (상호작용 p={pi:.3f})",
             color='#b71c1c' if 'RTM 아티팩트' in tag else '#1b5e20')
# 4th panel: interaction beta summary for all metrics
res=[(m.replace('d_',''),*inter_p(m)) for m in ['d_'+x for x in METS]]
res=sorted(res,key=lambda x:x[2])
ax=axes[3]; names=[r[0] for r in res]; ps=[r[2] for r in res]; betas=[r[1] for r in res]
cols=['#1b5e20' if p<0.05 else '#bdbdbd' for p in ps]
ax.barh(range(len(names)),betas,color=cols)
ax.set_yticks(range(len(names))); ax.set_yticklabels([f"{n} (p={p:.3f})" for n,p in zip(names,ps)],fontsize=9)
ax.axvline(0,color='k',ls='--'); ax.set_xlabel("impact×complexity 상호작용 beta")
ax.set_title("RTM 통제 후 상호작용\n(초록=유의 p<0.05)")
fig.suptitle("J05. RTM 통제: 초기 복잡도의 교란 완충효과 — 대조군과 기울기가 갈라지면 진짜(Deep_Gap·Gini), 평행이면 RTM(Canopy_Ht)",fontsize=12)
fig.tight_layout(); save(fig,"J05_rtm_survivors.png")

# ===== J06: calculation walkthrough on a real plot =====
# pick an impact plot with clear Deep_Gap change and >=2 pre & >=2 post surveys
cand=None
for pid in impmap:
    g=py[py.plotID==pid]; sy=impmap[pid]
    b=g[g.year<sy]; a=g[g.year>=sy]
    if len(b)>=2 and len(a)>=2 and abs(a['Deep_Gap'].mean()-b['Deep_Gap'].mean())>0.05:
        cand=pid; break
g=py[py.plotID==cand].sort_values('year'); sy=impmap[cand]
b=g[g.year<sy]; a=g[g.year>=sy]; pre=b['Deep_Gap'].mean(); post=a['Deep_Gap'].mean()
fig,ax=plt.subplots(1,2,figsize=(17,6))
# left: worked example time series
ax[0].plot(g['year'],g['Deep_Gap'],'o-',color='#37474f',ms=7,lw=1.5,label='Deep_Gap 관측(연도별)')
ax[0].axvline(sy,color='#d32f2f',ls='--',lw=2,label=f'교란연도 {int(sy)}')
ax[0].hlines(pre,b['year'].min(),sy,color='#1565c0',lw=3,label=f'교란 전 평균 = {pre:.3f}')
ax[0].hlines(post,sy,a['year'].max(),color='#e65100',lw=3,label=f'교란 후 평균 = {post:.3f}')
ax[0].annotate('',xy=(a['year'].max(),post),xytext=(a['year'].max(),pre),
   arrowprops=dict(arrowstyle='<->',color='green',lw=2))
ax[0].text(a['year'].max()+0.1,(pre+post)/2,f'Δ = {post-pre:+.3f}',color='green',fontsize=12,fontweight='bold',va='center')
ax[0].set_xlabel('연도'); ax[0].set_ylabel('Deep_Gap'); ax[0].legend(fontsize=9)
ax[0].set_title(f'① 계산: plot {cand}\nΔDeep_Gap = (교란 후 평균) - (교란 전 평균)')
# right: RTM-control schematic (impact vs control slope)
m='d_Deep_Gap'
for gflag,c,nm in [(0,'#1565c0','대조군(미교란): RTM 기준선'),(1,'#d32f2f','교란군')]:
    d=ch[ch.impact==gflag].dropna(subset=[m,'pre_complex'])
    sl,ic,r,p,se=st.linregress(d['pre_complex'],d[m]); xx=np.linspace(d.pre_complex.min(),d.pre_complex.max(),20)
    ax[1].plot(xx,ic+sl*xx,'-',color=c,lw=3,label=f'{nm} 기울기 {sl:+.3f}')
    ax[1].scatter(d['pre_complex'],d[m],s=10,alpha=.2,color=c)
bi,pi=inter_p(m)
ax[1].axhline(0,color='gray',ls=':'); ax[1].set_xlabel('초기 복잡도 (pre rugosity)'); ax[1].set_ylabel('ΔDeep_Gap')
ax[1].legend(fontsize=10); ax[1].set_title(f'② RTM 통제: 교란 vs 대조 기울기 비교\n상호작용(차이) = {bi:+.3f}, p={pi:.3f} -> 갈라짐=진짜 효과')
fig.suptitle("J06. 계산 방법 — ① plot별 교란 전/후 평균차(Δ) -> ② 교란군·대조군 기울기 차이(=RTM 제거한 순수 교란효과)",fontsize=13)
fig.tight_layout(); save(fig,"J06_calculation.png")
print("DONE")
