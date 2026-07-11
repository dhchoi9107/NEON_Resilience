# P3 — 생산성: 위성 녹색도 '혹형'의 함정과 진짜 생산성–다양성 관계

> ⚠ **2026-07-11 대전환**: DHI "생산성–다양성 혹형"은 **생산성 현상이 아니라 위성 녹색도(NDVI) 아티팩트**로 판명. 독립 GPP(MODIS+PML) 검증으로 반전됨. 폴더명(_hump)은 이력상 유지.

## 주제 (수정)
위성 대리변수로 생산성–다양성 관계를 검정하되, **겉보기 혹형이 진짜 생산성 신호인지 독립 GPP로 검증**. 결과: **DHI(Sentinel 녹색도) 혹형은 캐노피 밀도/NDVI 포화 아티팩트**이고, 실제 생산성(GPP)으로는 **종-에너지(단조 증가)** 관계.

## 핵심 결과
1. **겉보기 혹형**: Hill q1·q2 ~ DHI누적 = ∩ (다항식 p<0.001), forest type 통제에도 생존(`L04` β=−0.287). → 처음엔 "진짜 생산성 혹형"으로 오판.
2. ★★ **GPP 삼각검증으로 혹형 반전** (`L05·L06`, `dhi_gpp_{validation,triangulation}.csv`):
   - **MODIS GPP(MOD17) β=+0.22, PML GPP β=+0.17 — 둘 다 단조증가(혹형 아님).** DHI만 ∩.
   - **두 GPP는 서로 r=0.85 일치**, **DHI–GPP는 r=0.46~0.49**(중간)뿐 → **DHI가 아웃라이어.**
   - 혹형 하강부(고DHI=저다양)는 **밀집 상록 캐노피 NDVI 포화** 산물 (ABBY: NDVI 최고·GPP 중간·다양성 최저).
3. **진짜 관계 = 종-에너지(species-energy)**: 독립 GPP로 보면 **생산성↑ → 다양성↑ 단조**.
4. **통일 방법론 교훈**: 임령혹형(생물지리)·DHI혹형(녹색도) **둘 다 겉보기 혹형이 엄밀검증서 무너짐** → 위성 겉보기 곡선은 forest type + **독립 GPP**로 반드시 검증.
> 미완: NEON 타워 GPP(gold standard) 검증은 AmeriFlux 계정 필요(대기).

## 그림 (figures/)
- **L06** 삼각검증 headline (DHI 혹형 vs MODIS·PML 단조) ★ / **L05** DHI vs MODIS GPP
- L03 겉보기 혹형(DHI) / L04 forest type 통제(DHI서 생존) / L01·L02 비선형 / G01·G02 DHI

## 결과 (results/)
**dhi_gpp_validation**(DHI-MODIS 상관), **dhi_gpp_triangulation**(프록시별 2차부호), dhi_hump_confound, obj23_dhi_coeff, nonlinear, gam

## 데이터
`../../data/plot_dhi_sentinel.csv`, `plot_modis_gpp.csv`, `plot_pml_gpp.csv`

## 관련 스크립트 (../../scripts/obj23/)
`60·61`(DHI GEE)·`120~122`(비선형)·`147`(forest type)·**`149`(MODIS GPP)·`150`(PML GPP 삼각검증)**
