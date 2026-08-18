# -*- coding: utf-8 -*-
"""
Archive raw NEON source data from the external drive to NAS Volume2.

Why not cp/robocopy: robocopy failed against this // mount (rc=16, 2026-07-28) and
plain `cp -n` would silently skip a file that was truncated by an interrupted run.
This copies file-by-file and treats a destination file as done only if its size
matches the source exactly, so the job is safely resumable: re-running verifies
and fills gaps instead of trusting existence.

NOTHING IS DELETED. Deletion is a separate, explicit decision after verification.

Usage:  python archive_raw_to_nas.py [--verify-only]
Out:    scripts_pipeline/_pipeline_state/archive_raw_nas.log  (progress + summary)
"""
import os, sys, time, shutil, traceback

SRC_ROOT = r"E:\neon_lidar"
DST_ROOT = r"\\10.10.170.55\Volume2\NEON_raw"
FOLDERS  = ["DP1.30003.001", "vegetation_indices"]     # raw LiDAR point clouds, then raw VI tiles
LOG      = r"C:\Users\star1\Documents\GitHub\NEON_Resilience\scripts_pipeline\_pipeline_state\archive_raw_nas.log"
VERIFY_ONLY = "--verify-only" in sys.argv

def log(msg):
    line = f"[{time.strftime('%m-%d %H:%M:%S')}] {msg}"
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n"); fh.flush()
    print(line, flush=True)

def human(n):
    for u in ("B","KB","MB","GB","TB"):
        if n < 1024: return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} PB"

log("=" * 70)
log(f"START {'VERIFY-ONLY' if VERIFY_ONLY else 'COPY'}  src={SRC_ROOT}  dst={DST_ROOT}")

grand = dict(files=0, copied=0, skipped=0, failed=0, bytes_copied=0, bytes_total=0)
t0 = time.time()

for folder in FOLDERS:
    src_dir = os.path.join(SRC_ROOT, folder)
    dst_dir = os.path.join(DST_ROOT, folder)
    if not os.path.isdir(src_dir):
        log(f"SKIP {folder}: source not found"); continue
    log(f"--- {folder}: scanning ...")
    items = []
    for root, _dirs, files in os.walk(src_dir):
        for fn in files:
            if fn.startswith("._"):        # macOS resource forks, not data
                continue
            s = os.path.join(root, fn)
            try: sz = os.path.getsize(s)
            except OSError: continue
            items.append((s, os.path.join(dst_dir, os.path.relpath(s, src_dir)), sz))
    tot_bytes = sum(i[2] for i in items)
    grand["files"] += len(items); grand["bytes_total"] += tot_bytes
    log(f"--- {folder}: {len(items):,} files, {human(tot_bytes)}")

    done_bytes = 0; last_report = time.time()
    for i, (s, d, sz) in enumerate(items, 1):
        try:
            if os.path.exists(d) and os.path.getsize(d) == sz:
                grand["skipped"] += 1; done_bytes += sz
            elif VERIFY_ONLY:
                grand["failed"] += 1
                log(f"  MISSING/SIZE-MISMATCH: {os.path.relpath(s, SRC_ROOT)}")
            else:
                os.makedirs(os.path.dirname(d), exist_ok=True)
                tmp = d + ".part"
                shutil.copy2(s, tmp)
                if os.path.getsize(tmp) != sz:
                    os.remove(tmp); raise IOError("size mismatch after copy")
                os.replace(tmp, d)
                grand["copied"] += 1; grand["bytes_copied"] += sz; done_bytes += sz
        except Exception as e:
            grand["failed"] += 1
            log(f"  FAIL {os.path.relpath(s, SRC_ROOT)}: {type(e).__name__} {e}")
        if time.time() - last_report > 300:          # progress every 5 min
            el = time.time() - t0
            rate = grand["bytes_copied"] / el if el else 0
            remain = tot_bytes - done_bytes
            eta = remain / rate / 3600 if rate > 0 else float("nan")
            log(f"  {folder}: {i:,}/{len(items):,} files | {human(done_bytes)}/{human(tot_bytes)} "
                f"| copied {human(grand['bytes_copied'])} @ {human(rate)}/s | ETA {eta:.1f} h")
            last_report = time.time()
    log(f"--- {folder}: done (copied {grand['copied']:,}, verified-existing {grand['skipped']:,}, failed {grand['failed']:,})")

el = time.time() - t0
log("=" * 70)
log(f"SUMMARY  files={grand['files']:,}  copied={grand['copied']:,}  already-ok={grand['skipped']:,}  failed={grand['failed']:,}")
log(f"         transferred {human(grand['bytes_copied'])} of {human(grand['bytes_total'])} in {el/3600:.2f} h")
log("STATUS: " + ("ALL VERIFIED — safe to consider deletion" if grand["failed"] == 0
                  else f"{grand['failed']} PROBLEM FILES — DO NOT DELETE, re-run to fix"))
log("Source untouched; no files were deleted.")
