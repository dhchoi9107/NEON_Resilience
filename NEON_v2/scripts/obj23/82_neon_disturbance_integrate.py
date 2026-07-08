"""
Integrate NEON-NATIVE disturbance (authoritative) and cross-validate vs remote products.
Sources:
  (1) DP1.10111 sim_eventData  -> plot-level fire/flood/wind/drought/harvest/insect + date
  (2) vst plantStatus          -> in-situ stem mortality/damage severity (continuous)
  (3) vst remarks              -> text-mentioned disturbances
  (4) remote (MTBS/Hansen/NBR) -> independent corroboration
Output: plot_disturbance_neon.csv  (NEON-native primary disturbance flag/type/year/severity)
"""
import re, numpy as np, pandas as pd
D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
ev=pd.read_csv(D+r"\plot_neon_events.csv",low_memory=False)
ev['year']=pd.to_numeric(ev['startDate'].astype(str).str[:4],errors='coerce')

# classify disturbance type from methodTypeChoice / eventType (keep veg-relevant disturbances)
def dtype(row):
    s=f"{row.get('methodTypeChoice','')} {row.get('eventType','')}".lower()
    if 'fire' in s or 'burn' in s: return 'fire'
    if 'wind' in s or 'storm' in s or 'hurric' in s: return 'wind'
    if 'flood' in s: return 'flood'
    if 'drought' in s: return 'drought'
    if 'clearcut' in s or 'thinning' in s or 'harvest' in s or 'removal-mow' in s: return 'harvest'
    if 'populationspike' in s or 'insect' in s or 'invasive' in s: return 'insect'
    if 'naturaldisturbance' in s: return 'natural'
    return None
ev['dtype']=ev.apply(dtype,axis=1)
evd=ev.dropna(subset=['dtype','year']).copy()
print("disturbance events (veg-relevant):",len(evd),"| types:",evd.dtype.value_counts().to_dict())

# parse locationID -> affected basePlot IDs
def plots_of(s):
    if not isinstance(s,str): return []
    return list(set(re.findall(r'([A-Z]{4}_\d+)\.basePlot',s)))
rows=[]
for _,r in evd.iterrows():
    pl=plots_of(r['locationID']); n=len(pl)
    for pid in pl:
        rows.append((pid,int(r['year']),r['dtype'],n))
ev_plot=pd.DataFrame(rows,columns=['plotID','year','dtype','n_plots_event']).drop_duplicates()
# LOCALIZED event = touches <=12 basePlots (not a site-wide blanket record)
ev_loc=ev_plot[ev_plot['n_plots_event']<=12]
print("event->plot rows:",len(ev_plot),"| any-event plots:",ev_plot.plotID.nunique(),
      "| LOCALIZED-event plots:",ev_loc.plotID.nunique())

# per plot: most recent LOCALIZED disturbance event
evp=ev_loc.sort_values('year').groupby('plotID').agg(
    neon_event=('dtype','size'),neon_event_year=('year','max'),
    neon_event_type=('dtype','last')).reset_index()
evp['neon_event']=1

# in-situ severity + remarks
ins=pd.read_csv(D+r"\plot_insitu_disturbance.csv")
rmk=pd.read_csv(D+r"\plot_remarks_disturbance.csv")
rmk_p=rmk.groupby('plotID')['year'].max().reset_index().rename(columns={'year':'remark_year'}); rmk_p['remark']=1

# remote robust
rob=pd.read_csv(D+r"\plot_disturbance_robust.csv")[['plotID','disturbed','dist_year','severity','d_fire','d_loss','d_spec']]
rob=rob.rename(columns={'disturbed':'remote_dist','dist_year':'remote_year','severity':'remote_sev'})

# merge all on plot
m=ins.merge(evp,on='plotID',how='outer').merge(rmk_p,on='plotID',how='outer').merge(rob,on='plotID',how='outer')
for c in ['neon_event','remark','insitu_dist','remote_dist']: m[c]=m[c].fillna(0).astype(int)
m['insitu_spike']=m['insitu_spike'].fillna(0)
# NEON-native disturbance (primary) = documented LOCALIZED event OR mortality SPIKE OR burn OR remark
m['neon_dist']=((m.neon_event==1)|(m.insitu_spike>0.30)|(m.insitu_burn.fillna(0)>0)|(m.remark==1)).astype(int)
# NEON disturbance year (prefer documented event, else in-situ peak, else remark)
m['neon_dist_year']=m['neon_event_year'].fillna(m['insitu_peak_year']).fillna(m['remark_year'])
# severity: in-situ mortality/damage max (continuous ground truth)
m['neon_severity']=m['insitu_mortdmg_max'].fillna(0)
fb=np.where(m.insitu_burn.fillna(0)>0,'fire',np.where(m.insitu_insect.fillna(0)>0,'insect','mortality'))
m['neon_dist_type']=m['neon_event_type'].fillna(pd.Series(fb,index=m.index))
# BACI-compatible columns (dist_year ONLY for disturbed plots)
m['disturbed']=m['neon_dist']
m['dist_year']=np.where(m['neon_dist']==1,m['neon_dist_year'],np.nan)
m['severity']=m['neon_severity']; m['recency']=2025-m['dist_year']
m.to_csv(D+r"\plot_disturbance_neon.csv",index=False)

print(f"\n=== NEON-native disturbance: {int(m.neon_dist.sum())}/{len(m)} plots ===")
print("  유형:",m[m.neon_dist==1]['neon_dist_type'].value_counts().to_dict())
# ===== cross-validation vs remote =====
v=m.dropna(subset=['neon_dist','remote_dist'])
both=((v.neon_dist==1)&(v.remote_dist==1)).sum(); neon_only=((v.neon_dist==1)&(v.remote_dist==0)).sum()
rem_only=((v.neon_dist==0)&(v.remote_dist==1)).sum(); neither=((v.neon_dist==0)&(v.remote_dist==0)).sum()
print(f"\n=== NEON vs 원격(MTBS/Hansen/NBR) 교차검증 ===")
print(f"  둘다 교란 {both} | NEON만 {neon_only} | 원격만 {rem_only} | 둘다 정상 {neither}")
print(f"  일치율: {(both+neither)/len(v)*100:.0f}%")
# 검증: NEON 문서화 화재 plot이 실제로 in-situ 고사·원격 화재가 높은가
print("\n=== 검증: NEON 이벤트 유형별 in-situ 고사/손상 & 원격 일치 ===")
mm=m[m.neon_event==1]
for t in ['fire','harvest','wind','flood','insect']:
    sub=mm[mm.neon_event_type==t]
    if len(sub)<3: continue
    print(f"  {t:8s} n={len(sub):3d} | in-situ mort/dmg 중앙 {sub.insitu_mortdmg_max.median():.2f} | "
          f"원격도 교란 {int((sub.remote_dist==1).sum())}/{len(sub)}")
