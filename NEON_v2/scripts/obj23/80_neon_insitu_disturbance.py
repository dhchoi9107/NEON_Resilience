"""
NEON IN-SITU disturbance from vst plantStatus (ground-truth stem mortality/damage) + remarks.
Per plot-year: fraction of stems dead / damaged / burned / insect-damaged.
This is direct field evidence of disturbance — better than remote NBR drops.
Output: plot_year_insitu.csv (per plot-year fractions) + plot_insitu_disturbance.csv (plot summary)
"""
import numpy as np, pandas as pd, re
VS=r"E:\neon_lidar\vegetation_structure"; D=r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data"
ai=pd.read_csv(VS+r"\vst_apparentindividual.csv",low_memory=False,
               usecols=['plotID','date','plantStatus','growthForm'])
ai['year']=pd.to_numeric(ai['date'].astype(str).str[:4],errors='coerce')
ai=ai.dropna(subset=['year','plantStatus']); ai['year']=ai['year'].astype(int)
ai=ai[ai['growthForm'].astype(str).str.contains('tree|sapling|shrub',case=False,na=False)]
s=ai['plantStatus'].astype(str)
ai['dead']=s.str.contains('dead|downed|lost, presumed|broken bole',case=False)&~s.str.contains('Live',case=False)
ai['damaged']=s.str.contains('damaged|broken bole',case=False)&s.str.contains('Live',case=False)
ai['burned']=s.str.contains('burned',case=False)
ai['insect']=s.str.contains('insect',case=False)
ai['disease']=s.str.contains('disease',case=False)

g=ai.groupby(['plotID','year'])
py=g.agg(n=('plantStatus','size'),frac_dead=('dead','mean'),frac_dmg=('damaged','mean'),
         frac_burn=('burned','mean'),frac_insect=('insect','mean'),frac_dis=('disease','mean')).reset_index()
py['frac_mort_dmg']=py['frac_dead']+py['frac_dmg']   # total mortality+damage intensity (severity proxy)
py.to_csv(D+r"\plot_year_insitu.csv",index=False)
print("plot-year in-situ rows:",len(py),"| plots:",py.plotID.nunique())
print("frac_dead 중앙 %.2f | frac_dmg %.2f | burned>0 plot-yr: %d | insect>0: %d"%(
    py.frac_dead.median(),py.frac_dmg.median(),(py.frac_burn>0).sum(),(py.frac_insect>0).sum()))

# plot summary: peak mortality/damage year = in-situ disturbance event
rows=[]
for pid,gg in py.groupby('plotID'):
    gg=gg.sort_values('year')
    # baseline = min mort over surveys; event = max; spike = max-min
    peak=gg.loc[gg['frac_mort_dmg'].idxmax()]
    burn=gg['frac_burn'].max(); ins=gg['frac_insect'].max()
    spike=gg['frac_mort_dmg'].max()-gg['frac_mort_dmg'].min()
    rows.append((pid,float(gg['frac_mort_dmg'].max()),int(peak['year']),float(spike),
                 float(burn),float(ins),float(gg['frac_dead'].max())))
ins=pd.DataFrame(rows,columns=['plotID','insitu_mortdmg_max','insitu_peak_year','insitu_spike',
                               'insitu_burn','insitu_insect','insitu_dead_max'])
# flag in-situ disturbance: high mortality/damage OR any burn OR mortality spike
ins['insitu_dist']=((ins.insitu_mortdmg_max>0.4)|(ins.insitu_burn>0)|(ins.insitu_spike>0.25)).astype(int)
ins.to_csv(D+r"\plot_insitu_disturbance.csv",index=False)
print("\nplots:",len(ins),"| in-situ disturbed:",ins.insitu_dist.sum(),
      "(burned:%d, high-mort:%d)"%((ins.insitu_burn>0).sum(),(ins.insitu_mortdmg_max>0.4).sum()))

# remarks with disturbance keywords -> plot list
ppy=pd.read_csv(VS+r"\vst_perplotperyear.csv",low_memory=False,usecols=['plotID','date','remarks'])
ppy['year']=pd.to_numeric(ppy['date'].astype(str).str[:4],errors='coerce')
kw=ppy.dropna(subset=['remarks'])
kw=kw[kw['remarks'].astype(str).str.contains('fire|burn|harvest|logg|insect|beetle|wind|storm|hurricane|flood|killed|disturb|cut down|cleared',case=False)]
kw[['plotID','year','remarks']].to_csv(D+r"\plot_remarks_disturbance.csv",index=False)
print("remarks 교란 언급 plot-year:",len(kw),"| 고유 plot:",kw.plotID.nunique())
