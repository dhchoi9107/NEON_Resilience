# #2 (method-robustness) Formal variation partitioning (vegan-varpart style), domain-conditional.
# Cross-checks the manual semi-partial R^2 reported in the main framework: partitions the
# remote-sensing R^2 (beyond NEON domain) into unique structure / spectral / dynamics blocks
# plus their shared fraction, using set-complement R^2 differences. If this agrees with the
# main semi-partial numbers, the block conclusions are not an artifact of one partitioning recipe.
import os, pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
D="NEON_v2/data/"
OUTDIR="NEON_v2/papers/UNIFIED/results/O0_framework"; os.makedirs(OUTDIR, exist_ok=True)
S4=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"]; P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"]
f=pd.read_csv(D+"FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()
def prep(resp):
    d=f[[resp,'domain']+S4+P+Dn].dropna().copy()
    y=np.log(d[resp]) if 'Hill' in resp else d[resp].astype(float); d['y']=(y-y.mean())/y.std()
    for c in S4+P+Dn: d[c]=(d[c]-d[c].mean())/d[c].std()
    return d
def varpart(d):
    A,B,Cc=S4,P,Dn                               # Cc (not C) — C() is patsy's categorical op
    o=lambda pr: smf.ols('y ~ C(domain)'+('' if not pr else ' + '+' + '.join(pr)),d).fit().rsquared
    dom=smf.ols('y ~ C(domain)',d).fit().rsquared
    abc=o(A+B+Cc)-dom                             # total RS R^2 beyond domain
    uA=abc-(o(B+Cc)-dom); uB=abc-(o(A+Cc)-dom); uC=abc-(o(A+B)-dom)   # unique = full - complement
    return dict(n=int(len(d)),total=abc,uStruct=uA,uSpec=uB,uDyn=uC,shared=abc-(uA+uB+uC))
rows=[]
for resp in ['Hill_q1','Hill_q2','LCBD_turnover_rare','LCBD_nestedness_rare']:
    vp=varpart(prep(resp)); vp['response']=resp; rows.append(vp)
    print(f"  {resp:22}: total={vp['total']:.3f} | 구조={vp['uStruct']:.3f} 분광={vp['uSpec']:.3f} "
          f"동역학={vp['uDyn']:.3f} shared={vp['shared']:.3f}")
out=pd.DataFrame(rows)[['response','n','total','uStruct','uSpec','uDyn','shared']].round(4)
out.to_csv(os.path.join(OUTDIR,"varpart_robust.csv"),index=False)
print(f"\nsaved -> {OUTDIR}/varpart_robust.csv")
