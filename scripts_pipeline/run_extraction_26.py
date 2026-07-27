"""
Autonomous extraction orchestrator for the 7 NEW sites (per-site gated).
Download runs sequentially in SITE order; a site's LiDAR is COMPLETE once the next
site has any marker (download moved on) OR the site has all expected year-markers
OR the disk has been idle (no growth) for 30 min (last site / stall).
Per site: wait -> compute_plot_fsd_1m. Then VI download -> VI compute.
Stop-on-error, skip-if-done, resumable. Log: scripts_pipeline/_pipeline_state/master_log.txt
"""
import os, sys, time, subprocess, glob, shutil
ROOT = r"C:\Users\star1\Documents\GitHub\NEON_Resilience"
PY = r"C:/Users/star1/anaconda3/python.exe"
MARK = "E:/neon_lidar"
ORDER = ["DELA", "LENO", "UKFS", "YELL", "BONA", "DEJU", "HEAL"]
EXPECT = {"DELA": 8, "LENO": 9, "UKFS": 8, "YELL": 6, "BONA": 7, "DEJU": 7, "HEAL": 7}
MSTATE = os.path.join(ROOT, "scripts_pipeline", "_pipeline_state")
os.makedirs(MSTATE, exist_ok=True)
LOG = os.path.join(MSTATE, "master_log.txt")

def log(m):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
    try: print(line, flush=True)
    except Exception: pass          # never crash on console cp949 encoding
    with open(LOG, "a", encoding="utf-8") as f: f.write(line + "\n")

def done(tag): return os.path.exists(os.path.join(MSTATE, f".done_{tag}"))
def mark(tag): open(os.path.join(MSTATE, f".done_{tag}"), "w").close()
def nmark(site): return len(glob.glob(os.path.join(MARK, f".done_{site}_*")))
def dfree(): return shutil.disk_usage("E:/")[2]

def run(cmd, tag, env_extra=None):
    if done(tag): log(f"SKIP {tag}"); return True
    log(f"START {tag}: {cmd}")
    env = dict(os.environ); env["PYTHONIOENCODING"] = "utf-8"; env.update(env_extra or {})
    r = subprocess.run(cmd, shell=True, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    tail = "\n".join((r.stdout or "").splitlines()[-8:])
    if r.returncode == 0: log(f"OK {tag}\n{tail}"); mark(tag); return True
    log(f"FAIL {tag} (rc={r.returncode})\n{tail}"); return False

def wait_site(site, idx):
    """Block until this site's LiDAR download is complete (robust signals)."""
    if done(f"dl_{site}"): return True
    later = ORDER[idx + 1:]
    log(f"WAIT download: {site} (expect {EXPECT[site]} years)")
    last_free, idle_since = dfree(), time.time()
    while True:
        n = nmark(site)
        if any(nmark(s) > 0 for s in later):
            log(f"  {site} complete (download advanced to a later site); markers={n}"); mark(f"dl_{site}"); return True
        if n >= EXPECT[site]:
            log(f"  {site} complete ({n}/{EXPECT[site]} markers)"); mark(f"dl_{site}"); return True
        f = dfree()
        if last_free - f > 1e8:  # >0.1GB downloaded since last check -> active
            last_free, idle_since = f, time.time()
        elif time.time() - idle_since > 1800:  # 30 min no disk growth -> done/stalled
            log(f"  disk idle 30min; {site} has {n}/{EXPECT[site]} -> proceeding"); mark(f"dl_{site}"); return (n > 0)
        time.sleep(120)

def main():
    log("===== ORCHESTRATOR v2 START (per-site gating) =====")
    # Stage 1: per-site  wait-download -> LiDAR compute
    for i, s in enumerate(ORDER):
        if not wait_site(s, i):
            log(f"HALT: {s} has no LiDAR data"); sys.exit(1)
        if not run(f'"{PY}" compute/compute_plot_fsd_1m.py --site {s}', f"fsd_{s}"):
            log(f"HALT at LiDAR compute {s}"); sys.exit(1)
    log("----- all LiDAR compute done -----")
    # Stage 2: VI download (token file present; correct savepath default)
    if not run(f'"{PY}" download/neon_veg_indices_batch.py download --product VI',
               "vi_download", {"SAVEPATH": "E:/neon_lidar/vegetation_indices", "SITES_ENV": ",".join(ORDER)}):
        log("HALT at VI download"); sys.exit(1)
    # Stage 3: VI compute -> plot_vi_neon_brdf rows
    if not run(f'"{PY}" NEON_v2/scripts/01_extract_vi_neon_brdf.py', "vi_compute", {"SITES_ENV": ",".join(ORDER)}):
        log("HALT at VI compute"); sys.exit(1)
    log("===== DOWNLOAD+COMPUTE COMPLETE — supervised finalization next =====")

if __name__ == "__main__":
    main()
