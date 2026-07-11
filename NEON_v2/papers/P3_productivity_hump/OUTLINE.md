# P3 아웃라인 — 위성 녹색도 '혹형'의 함정 (재프레이밍)

## 제목(안)
- (주) "A satellite greenness 'productivity–diversity hump' is an artifact: independent GPP reveals a monotonic species–energy relationship"
- (부/국문) 위성 녹색도(DHI)의 겉보기 생산성–다양성 혹형은 아티팩트다 — 독립 GPP로 본 종-에너지 관계

## 타깃 저널
Remote Sensing of Environment / Global Ecology and Biogeography / Methods in Ecology and Evolution (방법론적 cautionary)

## 초록 (1문단 초안)
Sentinel-2 **DHI(누적 녹색도)** 와 CONUS 19개 산림 종다양성은 **혹형(inverted-U)** 을 보였고, forest type 통제에도 살아남아 고전적 생산성–다양성 혹형처럼 보였다. 그러나 **두 독립 위성 GPP(MODIS MOD17, PML-V2)로 검증하자 혹형이 사라지고 단조 증가로 반전**되었다(2차 β: DHI −0.29 vs MODIS +0.22, PML +0.17). 두 GPP는 서로 r=0.85로 일치했으나 DHI와는 r=0.46–0.49에 그쳐, **DHI가 아웃라이어**였다. 혹형의 하강부(고DHI–저다양)는 **밀집 상록 캐노피의 NDVI 포화**(고녹색도이나 중간 생산성·저다양성)에서 비롯됐다. 즉 겉보기 혹형은 위성 녹색도가 캐노피 밀도를 생산성으로 오독한 아티팩트이며, 실제 생산성으로는 **종-에너지(단조) 관계**가 나타난다. 위성 대리변수로 생태 곡선을 추론할 때 **독립 생산성 검증의 필요성**을 보인다.

## 핵심 주장(기여)
1. **겉보기 productivity–diversity 혹형이 위성 녹색도 아티팩트**임을 독립 GPP 2종으로 증명(혹형→단조 반전).
2. **DHI ≠ 생산성**: NDVI 포화·캐노피 밀도/상록성이 DHI를 부풀림(GPP와 r~0.5뿐, GPP끼리는 0.85).
3. **진짜 관계 = 종-에너지 단조 증가**.
4. ★ **cautionary methodology**: 위성 겉보기 곡선(혹형 등)은 반드시 (a)forest type (b)**독립 GPP**로 검증. 프로젝트 내 **임령혹형(생물지리 아티팩트)** 과 쌍을 이루는 방법론 교훈.

## Figure 순서 (main)
1. **L06** — Hill q1 vs DHI(∩) vs MODIS(단조) vs PML(단조) 3패널 ★headline
2. **L05** — DHI vs MODIS GPP 직접 비교(혹형 유무)
3. 프록시 상관 매트릭스(DHI-MODIS-PML) + ABBY 등 아웃라이어 사이트 표(NDVI 포화 예시)
4. L04 — DHI 혹형이 forest type엔 생존(왜 처음에 속았나)
- Supp: L01·L02(다항식·GAM), G01·G02, `dhi_gpp_triangulation.csv`

## 구성
- Intro: 위성 생산성 대리변수와 생산성–다양성 논쟁 / 겉보기 곡선의 검증 필요
- Methods: DHI(Sentinel), **MODIS MOD17·PML-V2 GPP(GEE)**, forest type 통제, 다항식/GAM, 3프록시 비교
- Results: (1)DHI 혹형·forest type 생존 (2)GPP 검증서 혹형 반전 (3)DHI 아웃라이어·NDVI 포화 기전 (4)진짜 종-에너지 관계
- Discussion: 위성 녹색도의 함정, NDVI 포화, 종-에너지 해석, 위성 생태추론의 검증 원칙

## 검증 완결 (3중)
- 위성 DHI(녹색도) → 위성 GPP 2종(MODIS·PML) → **실측 NEON 타워 GPP(16사이트)** 모두 혹형 부정.
- 실측 타워와 상관: **DHI 0.35(약) vs MODIS 0.67·PML 0.69** = DHI가 생산성 아님을 gold standard로 확정.

## 한계 (정직)
- 타워 GPP 3사이트 결측(ORNL 타워없음, SOAP·TEAK GPP미처리) → 타워 검정 n=16(site수준, 저검정력이나 부호 명확).
- MODIS/PML GPP도 모델(500m) — 단 서로 및 타워와 일치(0.67~0.85)하고 DHI만 갈림이 핵심.
- between-site 스케일, Sentinel 2016+.
