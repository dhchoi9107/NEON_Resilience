# NEON 타워 GPP (AmeriFlux FLUXNET-1F) — 다운로드 지침

이 폴더에 AmeriFlux FLUXNET-1F 파일을 풀어두면 `scripts/obj23/151_neon_tower_gpp.py`가 자동 처리합니다.

## 절차 (약 10분)
1. **계정 생성**: https://ameriflux.lbl.gov → *Register* (이메일 인증 + Data Use Policy 동의). 무료.
2. **데이터 다운로드**: 로그인 → *Download Data* →
   - **Network / Affiliation = "NEON"** 으로 필터 (NEON 전 사이트 한 번에 선택)
   - 또는 산림 사이트 개별 선택 (US-x** 코드)
   - **Data Product = "AmeriFlux FLUXNET-1F"** (⭐ GPP 포함; "BASE"는 GPP 없음)
   - Download 요청 → 이메일로 zip 링크 도착
3. **압축 해제**: 받은 zip들을 **이 폴더(`data/ameriflux/`)** 에 풀기.
   - 각 사이트 zip 안의 **연간 파일** `AMF_US-xXX_FLUXNET-1F_*_YY_*.csv` 이 핵심 (GPP_NT_VUT_REF).
4. **실행**: `python scripts/obj23/151_neon_tower_gpp.py`
   → 타워 GPP 추출 → 위성 3프록시(DHI/MODIS/PML)와 상관 + 다양성~타워GPP 혹형 검정 → L07.

## 필요한 사이트 (19개 산림, US-x 코드)
ABBY=US-xAB, BART=US-xBR, BLAN=US-xBL, GRSM=US-xGR, HARV=US-xHA, JERC=US-xJE, MLBS=US-xML,
OSBS=US-xSB(?), RMNP=US-xRM, SCBI=US-xSC, SERC=US-xSE, SOAP=US-xSP(?), STEI=US-xST, TALL=US-xTA,
TEAK=US-xTE, TREE=US-xTR, UNDE=US-xUN, WREF=US-xWR
> (?) 표시는 코드 확인 필요 — 스크립트가 매핑 안 된 파일을 출력하니 그때 맞추면 됨.

## 주의
- 타워는 사이트당 1개 → **타워 GPP는 사이트 단위(n≤19)**. 다양성~타워GPP는 사이트 수준 검정.
- FLUXNET-1F가 없는 사이트는 BASE(NEE)만 있을 수 있음 → 그 경우 파티셔닝 필요(별도).
- 이 폴더의 원자료는 gitignore 권장(용량). 결과 CSV만 커밋.
