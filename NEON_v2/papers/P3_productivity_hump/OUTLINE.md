# P3 아웃라인 — 생산성: 위성 DHI 생산성–다양성 혹형

## 제목(안)
- (주) "Satellite dynamic habitat indices reveal a robust productivity–diversity hump across temperate forests"
- (부/국문) 위성 DHI로 검출한 견고한 생산성–다양성 혹형, 그리고 겉보기 혹형의 함정

## 타깃 저널
Global Ecology and Biogeography / Ecography / Remote Sensing of Environment

## 초록 (1문단 초안)
Sentinel-2 **Dynamic Habitat Indices(DHI)** 를 생산성 대리변수로, 19개 NEON 온대림의 종다양성과의 관계를 검정했다. Hill q1·q2는 DHI누적과 **혹형(inverted-U, 다항식 p<0.001, GAM EDF 3.2)** — 중간 생산성에서 다양성 최대 — 를 보였고, 선형 +0.12는 착시였다. 이 혹형은 **forest type·상록성 통제에도 견고**(p<1e-5)하고 evergreen·deciduous **biome 내부에서도** 나타나, 단순 생물지리 아티팩트가 아님을 확인했다. 대조적으로 **임령–다양성의 겉보기 혹형은 forest type 통제 시 소멸(p=0.98)** 하는 아티팩트였다. 위성 생산성이 고전적 생산성–다양성 혹형을 between-site 거시생태 스케일에서 포착함을 보이며, 겉보기 혹형 해석 시 biome 통제의 필요성을 제시한다.

## 핵심 주장(기여)
1. **위성 DHI로 생산성–다양성 혹형 검출** — 선형은 착시, 다항식/GAM 필수.
2. **혹형은 진짜**: forest type·상록성 통제 생존 + biome 내부 존재.
3. ★ **겉보기 혹형의 함정**: 임령혹형(아티팩트) vs DHI혹형(진짜)의 대비 = 방법론적 경고. 위성 혹형은 forest type 통제로 검증하라.
4. **스케일**: between-site 거시생태(PDR 고전 스케일); 생산성–기후–biome 얽힘은 이 스케일의 표준 한계.

## Figure 순서 (main)
1. **L03** — 생산성–다양성 혹형(GAM) ★headline
2. **L04** — 견고성: forest type별 곡선(biome 내부에도 혹형) + 임령혹형 대비
3. **G02** — DHI 산점도(응답별)
- Supp: L01·L02(다항식·GAM 전체), G01(DHI forest), `dhi_hump_confound.csv`(통제별 2차 p 표)

## 구성
- Intro: 생산성–다양성 논쟁, 위성 생산성(DHI), 겉보기 혹형의 교란 위험(임령 사례)
- Methods: DHI 계산(Sentinel-2 GEE 서버사이드), 다항식·GAM, **forest type/상록성 통제 검증**, 임령혹형 대조군
- Results: (1)혹형 검출 (2)견고성(forest type·biome 내부) (3)임령혹형과 대비(아티팩트)
- Discussion: 거시생태 스케일 해석, DHI를 생산성 대리로, 방법론적 경고(biome 통제)

## 한계 (정직)
- **between-site 스케일**(site 통제 시 소멸 = 관계 고유 스케일, within-site 생산성 폭 좁음)
- 생산성 = 기후·biome과 얽힘(모든 PDR 공통), Sentinel 2016+ 기간 제한, 유효 N 사이트 수준
