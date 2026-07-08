"""
Fetch NEON DP1.10111 'Site management and event reporting' via API (documents disturbance/
management events per plot: fire, harvest, insect, treatment, etc.). Small text product.
Output: plot_neon_events.csv
"""
import sys, io, urllib.request, pandas as pd
sys.path.insert(0,r"C:\Users\star1\Documents\GitHub\NEON_Resilience")
from site_config import api_get, NEON_API, SITES
DP="DP1.10111.001"

def get_months(site):
    d=api_get(f"{NEON_API}/products/{DP}")["data"]
    for s in d["siteCodes"]:
        if s["siteCode"]==site: return s.get("availableMonths",[])
    return []

frames=[]
for site in SITES:
    for mo in get_months(site):
        try:
            files=api_get(f"{NEON_API}/data/{DP}/{site}/{mo}")["data"]["files"]
        except Exception: continue
        # ONLY the event-data table (sim_eventData); republished monthly -> dedup later
        cand=[f for f in files if "sim_eventData" in f["name"] and "basic" in f["name"] and f["name"].endswith(".csv")]
        for f in cand:
            try:
                raw=urllib.request.urlopen(f["url"],timeout=60).read()
                df=pd.read_csv(io.BytesIO(raw),low_memory=False); df["__site"]=site; frames.append(df)
            except Exception: continue
    print(f"{site}: done",flush=True)

ev=pd.concat(frames,ignore_index=True,sort=False) if frames else pd.DataFrame()
if 'uid' in ev.columns: ev=ev.drop_duplicates('uid')
elif len(ev): ev=ev.drop_duplicates()
print("\ntotal UNIQUE event rows:",len(ev))
if len(ev):
    print("cols:",[c for c in ev.columns][:25])
    # disturbance-type columns
    for c in ev.columns:
        if any(k in c.lower() for k in ('type','disturb','manage','categor')):
            vc=ev[c].dropna().value_counts()
            if 0<len(vc)<40: print(f"\n[{c}]:\n{vc.head(20).to_string()}")
    ev.to_csv(r"C:\Users\star1\Documents\GitHub\NEON_Resilience\NEON_v2\data\plot_neon_events.csv",index=False)
    print("\nsaved plot_neon_events.csv")
