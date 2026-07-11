# P1 아웃라인 — 패턴: RS로 종다양성 예측 + 조절

## 제목(안)
- (주) "Structural and spectral remote sensing complementarily predict tree diversity, but disturbance, land use, and succession modulate the link"
- (부/국문) 구조·분광 원격탐사의 상보적 종다양성 예측과 그 관계의 맥락 조절

## 타깃 저널
Remote Sensing of Environment / Global Ecology and Biogeography / Methods in Ecology and Evolution

## 초록 (1문단 초안)
19개 NEON 온대림에서 LiDAR 구조 다양성과 분광 식생지수(NEON BRDF)로 종다양성(알파 Hill q1/q2, 베타 LCBD turnover·nestedness)을 다중스케일 예측했다. **분광 VI는 알파를, 구조는 알파·베타를** 담당해 상보적이었고(알파 최상위 SAVI/EVI +0.38, 베타는 구조 전용), 시간 변화(trend)는 구조에서만 베타 정보를 추가했다. 이 RS–종 관계는 **교란(희귀종 손실·recency), 경관 파편화(−0.35), 천이단계(임령)** 에 의해 조절되었다. 특히 겉보기 **임령–다양성 혹형은 forest type 통제 시 소멸하는 생물지리 아티팩트**로, 위성 기반 다양성 추론에서 biome 통제의 필요성을 보였다.

## 핵심 주장(기여)
1. **구조 vs 분광 상보성**: 알파=분광 VI 최강, 베타=구조 전용. (기존 "구조>분광"보다 정밀)
2. **시간 trend 증분**: 구조 생장궤적이 mean 못 담는 베타 정보 추가(분광 0).
3. **조절 3축**: 교란·경관·임령이 RS–종 관계를 바꿈(관계는 고정 아님).
4. **방법론**: 임령·(대조군 P3의 생산성) 겉보기 혹형 중 임령혹형은 아티팩트 → forest type 통제 필수.

## Figure 순서 (main)
1. **K01** — 종다양성 ~ RS(분광 개별 VI + 구조), 4응답 forest ★
2. **F04** — 분산분해(domain>site>plot, 다중스케일)
3. **M01** — 시간 trend 증분(구조만 베타)
4. **H03**(+I01) — 교란이 RS–다양성 관계·조성 조절(recency, BACI)
5. **K02·K03** — 경관 이질성 직접·조절
6. **N10** — 임령–다양성 혹형 = 생물지리 아티팩트(forest type 통제 소멸)
- Supp: F01·F02·F05~F09, H01·H02, I02·I03, G03~G05, N09

## 구성
- Intro: RS로 종다양성 읽기 가능한가? 구조 vs 분광, 알파 vs 베타 / 그 관계가 맥락에 따라 변하나
- Methods: NEON 19사이트, LiDAR 구조지표, NEON .002 BRDF VI, Hill·LCBD, 혼합모델(domain fixed+site random), 조절항, forest type
- Results: (1)예측 구조/분광·알파/베타 (2)조절 교란·경관·임령 (3)혹형 아티팩트 검증
- Discussion: 상보성 메커니즘, "언제 RS가 통하나", biogeographic caution, NASA Obj1·2 연결

## 한계
분광 Sentinel/BRDF 연도 제한, 임령 사이트단위(유효 N≈19), 조성 pooled(시간동태는 BACI로 보완)
