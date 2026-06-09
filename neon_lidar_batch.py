#!/usr/bin/env python3
"""
NEON 항공 LiDAR 점군(DP1.30003.001) 일괄 다운로드 파이프라인
대상: 육상 생태계 사이트 (사막/건조지/툰드라/수생 제외), 2022~2025

모드:
  python3 neon_lidar_batch.py estimate   # 가용 여부 + 총 용량만 산정 (다운로드 X)
  python3 neon_lidar_batch.py download    # 실제 다운로드 (resumable)

환경변수:
  SAVEPATH=/Volumes/SSD/neon   저장 위치 (기본 ~/neon_data)
  NEON_TOKEN=...               (선택) NEON API 토큰 — 속도제한 완화
"""
import os, sys, json, time, glob, urllib.request, urllib.error

DPID  = "DP1.30003.001"
YEARS = list(range(2013, 2027))   # 2013~2026
# 연도 부분집합 병렬용: YEARS_ENV="2013,2024" 로 특정 연도만 처리
if os.environ.get("YEARS_ENV"):
    YEARS = [int(y.strip()) for y in os.environ["YEARS_ENV"].split(",") if y.strip()]
SAVEPATH = os.environ.get("SAVEPATH", os.path.expanduser("~/neon_data"))
# 토큰: 환경변수 우선 → 스크립트 옆 .neon_token → ~/.neon_token 순으로 읽음
TOKEN    = os.environ.get("NEON_TOKEN")
if not TOKEN:
    for _tf in (os.path.join(os.path.dirname(os.path.abspath(__file__)), ".neon_token"),
                os.path.expanduser("~/.neon_token")):
        if os.path.exists(_tf):
            TOKEN = open(_tf).read().strip()
            break

# ── 사이트 그룹 (생태계 분류) ────────────────────────────────
CONFIRMED = [  # 산림·초원·습지·사바나 (포함 확정)
    "ABBY","BART","BLAN","BONA","CLBJ","DCFS","DEJU","DELA","DSNY","GRSM",
    "HARV","JERC","KONA","KONZ","LENO","MLBS","NOGP","ORNL","OSBS","PUUM",
    "RMNP","SCBI","SERC","SOAP","STEI","TALL","TEAK","TREE","UKFS","UNDE",
    "WOOD","WREF","YELL",
]
DRY_FOREST   = ["GUAN","LAJA","SJER"]   # 건조 산림/사바나 → 포함
SEMIARID_GRS = ["CPER","OAES","STER"]   # 반건조 초원      → 제외(기본)
TRANSITION   = ["HEAL"]                 # 툰드라 전이대    → 제외(기본)

# 제외(참고용): 수생/사막/툰드라
AQUATIC = ["ARIK","BLUE","CUPE","GUIL","HOPB","LIRO","MCDI","MCRA","PRIN","REDB","SYCA","WLOU"]
DESERT  = ["JORN","MOAB","ONAQ","SRER"]
TUNDRA  = ["BARR","TOOL","NIWO"]

# ── 최종 대상 ────────────────────────────────────────────────
# 온대 산림 큐레이션 세트 (HARV/BART/BLAN/SCBI/WREF 유사). ORNL은 2022-25 비행 없어 제외, GRSM 대체.
SITES = sorted([
    # 온대 활엽·혼효 + PNW 침엽 (1차)
    "HARV","BART","UNDE","BLAN","SCBI","SERC","MLBS","ORNL","GRSM","WREF","ABBY",
    # 추가(2차): 서부 산악 침엽 + 남동부 소나무 + 북부 혼효 (YELL 제외)
    "SOAP","TEAK","RMNP","JERC","OSBS","TALL","STEI","TREE",
])
# 전체 육상(36개)으로 되돌리려면:
# SITES = sorted(CONFIRMED + DRY_FOREST)
# 병렬 실행: SITES_ENV="HARV,ABBY,..." 로 사이트 부분집합 지정
if os.environ.get("SITES_ENV"):
    SITES = [s.strip() for s in os.environ["SITES_ENV"].split(",") if s.strip()]

API = "https://data.neonscience.org/api/v0"

def _get(url):
    req = urllib.request.Request(url)
    if TOKEN:
        req.add_header("X-API-Token", TOKEN)
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

_PROD_CACHE = {}
def available_months(site):
    """해당 사이트의 점군 가용 'YYYY-MM' 목록 (release+provisional)."""
    if not _PROD_CACHE:
        d = _get(f"{API}/products/{DPID}")["data"]
        for s in d["siteCodes"]:
            _PROD_CACHE[s["siteCode"]] = sorted(s.get("availableMonths", []))
    return _PROD_CACHE.get(site, [])

def month_size_bytes(site, ym):
    """site/YYYY-MM 의 .laz 파일 총 바이트."""
    try:
        d = _get(f"{API}/data/{DPID}/{site}/{ym}")["data"]
    except urllib.error.HTTPError:
        return 0, 0
    # 분류 타일(.laz)만 집계, 미분류 flightline(Laz 폴더) 제외
    files = [f for f in d.get("files", [])
             if f["name"].lower().endswith(".laz")
             and "unclassified" not in f["name"].lower()]
    return sum(f.get("size", 0) for f in files), len(files)

def target_months(site):
    """2022~2025 에 해당하는 가용 월만."""
    return [ym for ym in available_months(site) if int(ym[:4]) in YEARS]

# ── ESTIMATE ─────────────────────────────────────────────────
def estimate():
    print(f"대상 사이트 {len(SITES)}개 × {YEARS[0]}~{YEARS[-1]}  용량 산정\n")
    print(f"{'site':5} {'months(2022-25)':22} {'tiles':>6} {'size':>10}")
    grand_b = grand_t = 0
    rows = []
    for site in SITES:
        yms = target_months(site)
        sb = st = 0
        for ym in yms:
            b, n = month_size_bytes(site, ym)
            sb += b; st += n
        grand_b += sb; grand_t += st
        rows.append((site, yms, st, sb))
        print(f"{site:5} {','.join(y[2:] for y in yms) or '-':22} {st:>6} {sb/1e9:>8.1f}G")
        time.sleep(0.1)
    print("-"*50)
    print(f"합계: 타일 {grand_t:,}개,  총 {grand_b/1e12:.2f} TB ({grand_b/1e9:.0f} GB)")
    print(f"\n저장 위치: {SAVEPATH}")
    return grand_b

# ── DOWNLOAD (resumable) ─────────────────────────────────────
def download():
    import neonutilities as nu
    # ── Laz(미분류 flightline) 제외 패치 ──────────────────────────
    # get_file_urls 결과에서 .../DiscreteLidar/Laz/ 파일을 빼버려
    # ClassifiedPointCloud(분류·타일) + Metadata 만 받음 → 용량/시간 ~49% 절약.
    # 나머지 로직(체크섬 재개·provisional 처리)은 그대로 유지.
    if os.environ.get("EXCLUDE_LAZ", "1") == "1":
        from neonutilities import aop_download as _aop
        _orig_gfu = _aop.get_file_urls
        def _no_laz(urls, token=None):
            df, rel = _orig_gfu(urls, token=token)
            if len(df):
                df = df[~df["url"].str.contains("/Laz/", na=False)].reset_index(drop=True)
            return df, rel
        _aop.get_file_urls = _no_laz
    os.makedirs(SAVEPATH, exist_ok=True)
    log = open(os.path.join(SAVEPATH, "download_log.txt"), "a")
    def L(m):
        print(m); log.write(m+"\n"); log.flush()
    for site in SITES:
        yrs = sorted({int(ym[:4]) for ym in target_months(site)})
        for yr in yrs:
            marker = os.path.join(SAVEPATH, f".done_{site}_{yr}")
            if os.path.exists(marker):
                L(f"[skip] {site} {yr} (이미 완료)"); continue
            L(f"[get ] {site} {yr} 다운로드 시작 {time.strftime('%H:%M:%S')}")
            try:
                nu.by_file_aop(
                    dpid=DPID, site=site, year=str(yr),
                    include_provisional=True,   # 2024~2025 미확정 포함
                    check_size=False,           # 일괄: 용량 프롬프트 비활성
                    savepath=SAVEPATH, token=TOKEN,
                    skip_if_exists=True,        # 체크섬 일치 파일만 건너뜀(재개)
                    overwrite="yes",            # 체크섬 불일치(잘린/변경) 파일은 자동 재다운로드 → 무결성 보장
                )
                open(marker, "w").close()
                L(f"[done] {site} {yr}")
            except Exception as e:
                L(f"[ERR ] {site} {yr}: {type(e).__name__}: {e}")
    log.close()
    verify()

# ── VERIFY ───────────────────────────────────────────────────
def verify():
    import laspy
    laz = glob.glob(os.path.join(SAVEPATH, "**", "*.laz"), recursive=True)
    print(f"\n검증: .laz {len(laz)}개")
    if not laz:
        return
    f0 = sorted(laz)[0]
    with laspy.open(f0) as fh:
        h = fh.header
        print(f"  예시 {os.path.basename(f0)} | pts={h.point_count:,} | CRS={h.parse_crs()}")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "estimate"
    {"estimate": estimate, "download": download, "verify": verify}[mode]()
