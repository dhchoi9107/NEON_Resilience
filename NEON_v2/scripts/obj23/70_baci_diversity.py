"""
BACI step 1 — time-resolved taxonomic diversity, BEFORE vs AFTER a split year, per plot.
Impact plots split at their disturbance year; control plots at their site's median impact year.
Pools live woody individuals within each period -> effort-robust Hill q1/q2 + sample coverage.
Output: plot_baci_diversity.csv  (plotID, impact, severity, recency,
        rich_before/after, q1_before/after, q2_before/after, cov_before/after, dq1, dq2, drich)
"""
import numpy as np, pandas as pd, sys
VS=r"E:\neon_lidar\vegetation_structure"
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
DISTFILE=sys.argv[1] if len(sys.argv)>1 else "plot_disturbance_robust.csv"
OUTSUF=sys.argv[2] if len(sys.argv)>2 else ""

ai=pd.read_csv(VS+r"\vst_apparentindividual.csv",low_memory=False,
               usecols=['plotID','individualID','date','plantStatus','growthForm'])
mt=pd.read_csv(VS+r"\vst_mappingandtagging.csv",low_memory=False,
               usecols=['individualID','taxonID','scientificName'])
ai['year']=pd.to_numeric(ai['date'].astype(str).str[:4],errors='coerce')
ai=ai[ai['plantStatus'].astype(str).str.contains('Live',case=False,na=False)]
ai=ai[ai['growthForm'].astype(str).str.contains('tree|sapling',case=False,na=False)]
sp=ai.merge(mt.drop_duplicates('individualID'),on='individualID',how='left').dropna(subset=['taxonID','year'])
sp['year']=sp['year'].astype(int)
print("live woody stems w/ species & year:",len(sp),"| plots:",sp.plotID.nunique())

rob=pd.read_csv(D+"\\"+DISTFILE)
imp=rob[rob.dist_year.notna()][['plotID','dist_year','severity','recency']]
# site of each plot
site=sp.groupby('plotID').first().reset_index()[['plotID']]
sp['siteID']=sp['plotID'].str[:4]
# split year: impact -> own dist_year; control -> site median impact year (fallback global)
impmap=imp.set_index('plotID')['dist_year'].to_dict()
site_med=imp.assign(siteID=imp.plotID.str[:4]).groupby('siteID')['dist_year'].median().to_dict()
gmed=float(imp['dist_year'].median())
def split_year(pid):
    if pid in impmap: return impmap[pid]
    return site_med.get(pid[:4],gmed)
sp['split']=sp['plotID'].map(split_year)
sp['period']=np.where(sp['year']<sp['split'],'before','after')

def hill_cov(counts):
    c=np.asarray(counts); c=c[c>0]; N=c.sum()
    if N<3 or len(c)<1: return np.nan,np.nan,np.nan,np.nan
    p=c/N
    q1=float(np.exp(-(p*np.log(p)).sum()))          # exp(Shannon)
    q2=float(1.0/(p**2).sum())                       # inverse Simpson
    f1=int((c==1).sum()); f2=int((c==2).sum())
    cov=1-(f1/N)*(((N-1)*f1)/((N-1)*f1+2*f2)) if f1>0 and ((N-1)*f1+2*f2)>0 else 1-0/N
    return float(len(c)),q1,q2,float(cov)

rows=[]
for (pid,per),g in sp.groupby(['plotID','period']):
    counts=g.groupby('taxonID').size().values
    rich,q1,q2,cov=hill_cov(counts)
    rows.append((pid,per,rich,q1,q2,cov,g['year'].nunique()))
pp=pd.DataFrame(rows,columns=['plotID','period','rich','q1','q2','cov','nyr'])
wide=pp.pivot(index='plotID',columns='period',values=['rich','q1','q2','cov','nyr'])
wide.columns=[f"{a}_{b}" for a,b in wide.columns]; wide=wide.reset_index()
# require both periods present
both=wide.dropna(subset=['q1_before','q1_after']).copy()
for m in ['rich','q1','q2']:
    both[f'd{m}']=both[f'{m}_after']-both[f'{m}_before']
both['impact']=both['plotID'].isin(impmap).astype(int)
both=both.merge(imp[['plotID','severity','recency']],on='plotID',how='left')
both['severity']=both['severity'].fillna(0);
both.to_csv(D+"\\plot_baci_diversity"+OUTSUF+".csv",index=False)
print(f"\nplots with before&after diversity: {len(both)} | impact {both.impact.sum()} / control {(both.impact==0).sum()}")
print("ΔHillq1: impact median %.2f | control median %.2f"%(
    both[both.impact==1]['dq1'].median(), both[both.impact==0]['dq1'].median()))
print("coverage before/after median: %.2f / %.2f"%(both['cov_before'].median(),both['cov_after'].median()))
