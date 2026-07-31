# #1 Bayesian multilevel — divergence-resolved: target_accept=0.99, more tuning,
# HalfNormal priors on group sigmas (site var ~0 causes the funnel). Report divergences + rhat.
# Saves a full coefficient table (all predictors, 95% HDI, prob-of-direction) + a meta row per
# response (divergences, max rhat, group sigmas) to results/O0_framework/ so the manuscript can
# cite the Bayesian inference as the primary test for dynamics -> nestedness.
import os, pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
import bambi as bmb, arviz as az
D="NEON_v2/data/"
OUTDIR="NEON_v2/papers/UNIFIED/results/O0_framework"; os.makedirs(OUTDIR, exist_ok=True)
f=pd.read_csv(D+"FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];keep=S+P+Dn
def block(c): return "structure" if c in S else "spectral" if c in P else "dynamics"
coef_rows=[]; meta_rows=[]
for resp in ['Hill_q1','LCBD_nestedness_rare']:
    d=f[[resp,'domain','siteID']+keep].dropna().copy()
    y=np.log(d[resp]) if 'Hill' in resp else d[resp].astype(float); d['y']=(y-y.mean())/y.std()
    for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
    priors={"1|domain":bmb.Prior("Normal",mu=0,sigma=bmb.Prior("HalfNormal",sigma=0.5)),
            "1|domain:siteID":bmb.Prior("Normal",mu=0,sigma=bmb.Prior("HalfNormal",sigma=0.3))}
    m=bmb.Model('y ~ '+' + '.join(keep)+' + (1|domain) + (1|domain:siteID)',d,priors=priors)
    idata=m.fit(draws=2000,tune=3000,chains=4,cores=1,target_accept=0.99,progressbar=False,random_seed=11)
    ndiv=int(idata.sample_stats.diverging.values.sum()); post=idata.posterior
    rh=az.rhat(idata,var_names=keep); maxrh=float(max(float(rh[c].values) for c in keep if c in rh))
    # group sigmas (funnel diagnostic: site sigma ~0)
    def gsig(name):
        vn=[v for v in post.data_vars if name in v and "sigma" in v]
        return float(post[vn[0]].values.mean()) if vn else np.nan
    sig_dom, sig_site = gsig("1|domain"), gsig("1|domain:siteID")
    print(f"\n===== {resp}: divergences={ndiv} (was 117) | max rhat={maxrh:.3f} | sigma(domain)={sig_dom:.3f} sigma(site)={sig_site:.3f} =====")
    meta_rows.append({"response":resp,"n":int(len(d)),"divergences":ndiv,"max_rhat":round(maxrh,4),
                      "sigma_domain":round(sig_dom,4),"sigma_site":round(sig_site,4)})
    s_ex=d_ex=0
    for c in keep:
        v=post[c].values.flatten(); lo,hi=np.percentile(v,2.5),np.percentile(v,97.5); ex=bool(lo>0 or hi<0)
        pd_dir=float(max((v>0).mean(),(v<0).mean()))   # posterior probability of direction
        if c in S and ex: s_ex+=1
        if c in Dn and ex: d_ex+=1
        coef_rows.append({"response":resp,"predictor":c,"block":block(c),"beta":round(float(v.mean()),4),
                          "hdi_2.5":round(float(lo),4),"hdi_97.5":round(float(hi),4),
                          "pd":round(pd_dir,4),"credible":ex})
        if ex: print(f"  {c:16} {v.mean():+.3f} [{lo:+.3f},{hi:+.3f}]  pd={pd_dir:.3f}")
    print(f"  블록: 구조 {s_ex}/4 0제외 | 동역학 {d_ex}/7 0제외")
pd.DataFrame(coef_rows).to_csv(os.path.join(OUTDIR,"bayes_multilevel_coeffs.csv"),index=False)
pd.DataFrame(meta_rows).to_csv(os.path.join(OUTDIR,"bayes_multilevel_meta.csv"),index=False)
print(f"\nsaved -> {OUTDIR}/bayes_multilevel_coeffs.csv & _meta.csv")
