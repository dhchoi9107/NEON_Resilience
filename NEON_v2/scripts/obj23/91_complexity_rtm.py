"""
Rigorous test of Choi-2023 'complexity offset' controlling for REGRESSION TO THE MEAN (RTM).

Two RTM controls:
 (A) CONTROL PLOTS as RTM baseline. RTM affects impact & control equally, so the real
     disturbance-specific complexity effect = impact x complexity INTERACTION
     (Δ ~ impact * z(complexity) + severity + (1|site)). The complexity MAIN effect is RTM.
 (B) Mathematical-coupling fix: regress Δ on the (pre+post)/2 MEAN complexity instead of pre
     (regressing Δ=post-pre on pre alone induces a spurious ~-0.5 correlation).

If the impact:complexity interaction is ~0 / n.s. -> the complexity-Δ relationship is just RTM,
NOT genuine buffering. If significant in buffering direction -> real complexity offset.
Output: results/complexity_rtm.csv + figure J04
"""
import os, numpy as np, pandas as pd, scipy.stats as st
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2"
D,F,R=os.path.join(ROOT,"data"),os.path.join(ROOT,"figures"),os.path.join(ROOT,"results")
def z(s): return (s-s.mean())/s.std()
def save(fig,n): fig.savefig(os.path.join(F,n),dpi=118,bbox_inches='tight'); plt.close(fig); print(" saved",n)

neon=pd.read_csv(os.path.join(D,"plot_disturbance_neon.csv"))
imp=neon[neon.disturbed==1].dropna(subset=['dist_year'])
impmap=imp.set_index('plotID')['dist_year'].to_dict()
sev=neon.set_index('plotID')['severity'].to_dict()
imp['siteID']=imp.plotID.str[:4]; smed=imp.groupby('siteID')['dist_year'].median().to_dict()
gmed=float(imp['dist_year'].median())
def split(pid): return impmap.get(pid, smed.get(pid[:4],gmed))

py=pd.read_csv(os.path.join(D,"per_year_v2.csv"))
MET=['Canopy_Ht','FHD','Rugosity']
rows=[]
for pid,g in py.groupby('plotID'):
    sy=split(pid); b=g[g.year<sy]; a=g[g.year>=sy]
    if len(b)<1 or len(a)<1: continue
    rec={'plotID':pid,'impact':int(pid in impmap),'severity':sev.get(pid,0)}
    for m in MET:
        if m in g:
            pre=b[m].mean(); post=a[m].mean()
            rec['d_'+m]=post-pre; rec['pre_'+m]=pre; rec['mean_'+m]=(pre+post)/2
    rows.append(rec)
ch=pd.DataFrame(rows); ch['siteID']=ch.plotID.str[:4]
ch['severity']=ch['severity'].fillna(0)
bac=pd.read_csv(os.path.join(D,"plot_baci_diversity_neon.csv"))[['plotID','drich','dq1']]
ch=ch.merge(bac,on='plotID',how='left')
print(f"plots: {len(ch)} | impact {int(ch.impact.sum())} control {int((ch.impact==0).sum())}")

res=[]
def lmm(d,form,term):
    d=d.dropna(subset=[c for c in form.replace('~',' ').replace('*',' ').replace('+',' ').split() if c in d]+[form.split('~')[0].strip()])
    try:
        m=smf.mixedlm(form,d,groups=d['siteID']).fit(reml=True,method='lbfgs')
        return m.fe_params.get(term,np.nan),m.pvalues.get(term,np.nan),len(d)
    except Exception: return np.nan,np.nan,len(d)

print("\n=== ΔCanopy_Ht / ΔFHD / ΔRichness: naive vs RTM-controlled ===")
COMPLEX='Rugosity'  # initial complexity metric
for dv,lab in [('d_Canopy_Ht','ΔCanopyHt'),('d_FHD','ΔFHD'),('drich','ΔRichness')]:
    d=ch.dropna(subset=[dv,'pre_'+COMPLEX]).copy()
    d['zpre']=z(d['pre_'+COMPLEX]); d['zmean']=z(d['mean_'+COMPLEX])
    # naive (impact plots only, pre): the biased estimate
    di=d[d.impact==1]; rn,pn=st.pearsonr(di['zpre'],di[dv]) if len(di)>10 else (np.nan,np.nan)
    # RTM baseline in controls
    dc=d[d.impact==0]; rc,pc=st.pearsonr(dc['zpre'],dc[dv]) if len(dc)>10 else (np.nan,np.nan)
    # (A) interaction: impact x complexity (real buffering net of RTM)
    bi,pi,n=lmm(d,f"{dv} ~ impact * zpre + severity",'impact:zpre')
    # (B) coupling-free: interaction using MEAN complexity
    bm,pm,_=lmm(d,f"{dv} ~ impact * zmean + severity",'impact:zmean')
    print(f"\n{lab}:")
    print(f"  naive(impact,pre): r={rn:+.2f} (p={pn:.3f})   <- biased by RTM+coupling")
    print(f"  RTM baseline(control,pre): r={rc:+.2f} (p={pc:.3f})   <- pure RTM slope")
    print(f"  (A) impact×complexity interaction: beta={bi:+.3f} (p={pi:.3f})  <- REAL buffering net of RTM")
    print(f"  (B) impact×MEANcomplexity (coupling-free): beta={bm:+.3f} (p={pm:.3f})")
    res.append(dict(response=lab,naive_r=rn,control_rtm_r=rc,inter_pre_beta=bi,inter_pre_p=pi,
                    inter_mean_beta=bm,inter_mean_p=pm,n=n))
pd.DataFrame(res).to_csv(os.path.join(R,"complexity_rtm.csv"),index=False)

# ---- J04: Δ vs complexity, impact vs control slopes ----
fig,axes=plt.subplots(1,3,figsize=(16,5))
for ax,(dv,lab) in zip(axes,[('d_Canopy_Ht','ΔCanopy height'),('d_FHD','ΔFHD'),('drich','ΔRichness')]):
    for g,c,nm in [(0,'#1565c0','control (RTM baseline)'),(1,'#d32f2f','impact')]:
        d=ch[ch.impact==g].dropna(subset=[dv,'pre_'+COMPLEX])
        if len(d)<8: continue
        ax.scatter(d['pre_'+COMPLEX],d[dv],s=14,alpha=.35,color=c)
        sl,ic,r,p,se=st.linregress(d['pre_'+COMPLEX],d[dv]);xx=np.linspace(d['pre_'+COMPLEX].min(),d['pre_'+COMPLEX].max(),20)
        ax.plot(xx,ic+sl*xx,'-',color=c,lw=2.3,label=f"{nm}: slope={sl:+.2f}")
    ax.axhline(0,color='gray',ls=':');ax.set_xlabel("initial complexity (pre rugosity)");ax.set_ylabel(lab)
    ax.legend(fontsize=8);ax.set_title(lab)
fig.suptitle("J04. RTM control — if impact slope ≈ control slope, complexity effect is just RTM (not buffering)",fontsize=12)
fig.tight_layout(); save(fig,"J04_rtm_control.png")
print("\nDONE")
