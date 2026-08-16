# Prior-sensitivity check for the Bayesian multilevel models (Supplement):
# group-scale HalfNormal priors doubled (domain 0.5->1.0, site 0.3->0.6).
# If the credible set and coefficient means are stable, the headline
# dynamics->nestedness result is not prior-driven.
import os, pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
import bambi as bmb, arviz as az
D="NEON_v2/data/"; OUT="NEON_v2/papers/UNIFIED/results/O0_framework"
f=pd.read_csv(D+"FINAL_v2_pooled_26.csv"); f=f[f.sample_coverage>=0.9].copy()
S=["Rugosity_mean","Vert_CV_mean","VCI_mean","LAI_mean"];P=["EVI_mean"]
Dn=["Rumple_trend","Vert_SD_trend","Vert_CV_trend","VCI_trend","FHD_trend","LAI_trend","Ht_Ratio_trend"];keep=S+P+Dn
rows=[]
for resp in ['Hill_q1','LCBD_nestedness_rare']:
    d=f[[resp,'domain','siteID']+keep].dropna().copy()
    y=np.log(d[resp]) if 'Hill' in resp else d[resp].astype(float); d['y']=(y-y.mean())/y.std()
    for c in keep: d[c]=(d[c]-d[c].mean())/d[c].std()
    priors={"1|domain":bmb.Prior("Normal",mu=0,sigma=bmb.Prior("HalfNormal",sigma=1.0)),
            "1|domain:siteID":bmb.Prior("Normal",mu=0,sigma=bmb.Prior("HalfNormal",sigma=0.6))}
    m=bmb.Model('y ~ '+' + '.join(keep)+' + (1|domain) + (1|domain:siteID)',d,priors=priors)
    idata=m.fit(draws=2000,tune=3000,chains=4,cores=1,target_accept=0.99,progressbar=False,random_seed=11)
    ndiv=int(idata.sample_stats.diverging.values.sum()); post=idata.posterior
    print(f"===== {resp} (2x priors): divergences={ndiv} =====",flush=True)
    for c in keep:
        v=post[c].values.flatten(); lo,hi=np.percentile(v,2.5),np.percentile(v,97.5)
        rows.append({"response":resp,"predictor":c,"beta_2x":round(float(v.mean()),4),
                     "hdi_2.5_2x":round(float(lo),4),"hdi_97.5_2x":round(float(hi),4),
                     "credible_2x":bool(lo>0 or hi<0),"divergences":ndiv})
pd.DataFrame(rows).to_csv(os.path.join(OUT,"bayes_prior_sensitivity.csv"),index=False)
print("saved -> bayes_prior_sensitivity.csv")
