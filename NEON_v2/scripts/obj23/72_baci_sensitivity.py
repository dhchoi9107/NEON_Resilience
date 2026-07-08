"""
BACI sensitivity — is the disturbance richness-loss robust to design choices?
Varies: (1) control split-year rule, (2) min individuals/period, (3) disturbance definition.
Re-pools individuals -> Hill/richness per period -> BACI mixedlm ΔDiv ~ impact + C(domain).
Output: results/baci_sensitivity.csv
"""
import numpy as np, pandas as pd, statsmodels.formula.api as smf, scipy.stats as st
VS=r"E:\neon_lidar\vegetation_structure"; D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
DOM={"HARV":"D01","BART":"D01","SCBI":"D02","SERC":"D02","BLAN":"D02","GRSM":"D07","MLBS":"D07","ORNL":"D07",
 "JERC":"D03","OSBS":"D03","TALL":"D08","UNDE":"D05","STEI":"D05","TREE":"D05","WREF":"D16","ABBY":"D16",
 "SOAP":"D17","TEAK":"D17","RMNP":"D10"}
# species x plot x year (once)
ai=pd.read_csv(VS+r"\vst_apparentindividual.csv",low_memory=False,usecols=['plotID','individualID','date','plantStatus','growthForm'])
mt=pd.read_csv(VS+r"\vst_mappingandtagging.csv",low_memory=False,usecols=['individualID','taxonID'])
ai['year']=pd.to_numeric(ai['date'].astype(str).str[:4],errors='coerce')
ai=ai[ai['plantStatus'].astype(str).str.contains('Live',case=False,na=False)]
ai=ai[ai['growthForm'].astype(str).str.contains('tree|sapling',case=False,na=False)]
sp=ai.merge(mt.drop_duplicates('individualID'),on='individualID',how='left').dropna(subset=['taxonID','year'])
sp['year']=sp['year'].astype(int); sp['siteID']=sp['plotID'].str[:4]

rob=pd.read_csv(D+r"\plot_disturbance_robust.csv")
imp_all=rob[rob.dist_year.notna()].copy(); imp_all['siteID']=imp_all.plotID.str[:4]
auth=rob[(rob.d_fire.fillna(0)+rob.d_loss.fillna(0))>0]   # authoritative only

def hill(counts):
    c=np.asarray(counts); c=c[c>0]; N=c.sum()
    if N<1: return np.nan,np.nan
    p=c/N; return float(len(c)), float(np.exp(-(p*np.log(p)).sum()))

def run_spec(impset, split_rule, min_N, seed=0):
    impset=impset.copy(); impset['siteID']=impset['plotID'].str[:4]
    impmap=impset.set_index('plotID')['dist_year'].to_dict()
    smed=impset.groupby('siteID')['dist_year'].median().to_dict(); gmed=float(impset['dist_year'].median())
    rng=np.random.RandomState(seed); years=impset['dist_year'].values
    ctrl_fixed={}
    def split(pid):
        if pid in impmap: return impmap[pid]
        if split_rule=='site_median': return smed.get(pid[:4],gmed)
        if split_rule=='global_median': return gmed
        if split_rule=='fixed_2018': return 2018
        if split_rule=='random_impact':
            if pid not in ctrl_fixed: ctrl_fixed[pid]=float(rng.choice(years))
            return ctrl_fixed[pid]
    rows=[]
    for pid,g in sp.groupby('plotID'):
        sy=split(pid); b=g[g.year<sy]; a=g[g.year>=sy]
        if len(b)<min_N or len(a)<min_N: continue
        rb,q1b=hill(b.groupby('taxonID').size().values); ra,q1a=hill(a.groupby('taxonID').size().values)
        rows.append((pid,pid[:4],int(pid in impmap),ra-rb,q1a-q1b))
    df=pd.DataFrame(rows,columns=['plotID','siteID','impact','drich','dq1'])
    out={}
    for dv in ['drich','dq1']:
        d=df.dropna(subset=[dv])
        imp_v=d[d.impact==1][dv]; ctl_v=d[d.impact==0][dv]
        # proper model for nested structure: site random intercept (no domain fixed -> avoids singular)
        try:
            mm=smf.mixedlm(f"{dv} ~ impact",d,groups=d['siteID']).fit(reml=True,method='lbfgs')
            beta,p=mm.fe_params['impact'],mm.pvalues['impact']
        except Exception:
            beta,p=np.nan,np.nan
        mwp=st.mannwhitneyu(imp_v,ctl_v).pvalue if len(imp_v)>3 and len(ctl_v)>3 else np.nan
        out[dv]=(beta,p,imp_v.median()-ctl_v.median(),mwp)
    return int(df.impact.sum()), int((df.impact==0).sum()), out

specs=[("all, site_median, N>=3", imp_all,'site_median',3),
       ("all, global_median, N>=3", imp_all,'global_median',3),
       ("all, fixed_2018, N>=3", imp_all,'fixed_2018',3),
       ("all, random_impact, N>=3", imp_all,'random_impact',3),
       ("all, site_median, N>=5", imp_all,'site_median',5),
       ("all, site_median, N>=10", imp_all,'site_median',10),
       ("AUTH only, site_median, N>=3", auth,'site_median',3),
       ("AUTH only, global_median, N>=3", auth,'global_median',3)]
res=[]
print(f"{'spec':32s} | nImp/nCtrl | ΔRich: site-lmm beta(p) | Δmed | MannW p")
for name,iset,rule,mn in specs:
    ni,nc,o=run_spec(iset,rule,mn)
    rb,rp,rmed,rmw=o['drich']
    sig='*' if rp<0.05 else ('.' if rp<0.1 else ' ')
    print(f"{name:32s} | {ni:3d}/{nc:3d}   | {rb:+.2f} (p={rp:.3f}){sig}        | {rmed:+.1f} | {rmw:.3f}")
    res.append(dict(spec=name,n_impact=ni,n_control=nc,drich_lmm_beta=rb,drich_lmm_p=rp,
                    drich_med_diff=rmed,drich_mannw_p=rmw,dq1_lmm_beta=o['dq1'][0],dq1_lmm_p=o['dq1'][1]))
pd.DataFrame(res).to_csv(D+r"\..\results\baci_sensitivity.csv",index=False)
print("\n→ ΔRich 음수 일관 & MannW<0.05 다수면 robust. DONE")
