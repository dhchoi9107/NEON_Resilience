# Concept Note — Canopy structural development *rate* from repeat airborne LiDAR

*(2026-08 초안. 다양성 논문(UNIFIED)에서 분리한 독립 논문 후보. 1-page scope.)*

## 제목 후보
- **"Direct measurement of canopy structural development rates from continental repeat airborne LiDAR, and their scaling with stand age"**
- (대안) "How fast do forests build structure? Metric-specific developmental timescales from repeat NEON LiDAR"

## 한 줄 요지
같은 플롯을 **반복 항공 LiDAR로 여러 해 관측**해 캐노피 구조 발달 *속도*(연간 추세)를 **직접 측정**하고, 그 속도가 **임령에 따라 지표별로 어떻게 감쇠**하는지, 그리고 **GPP 추세(탄소 축적)와 어떻게 결합**되는지를 대륙 규모(26 NEON 사이트, 온대~boreal)에서 정량화한다.

## 갭 (왜 새로운가)
- 기존 산림 천이·구조발달 연구는 대부분 **공간대치(space-for-time, 크로노시퀀스)** — 나이 다른 숲들을 한 시점에 비교해 발달을 *추론*한다.
- 여기서는 **동일 플롯의 반복 측정**으로 발달 속도를 *직접 관측*한다. 크로노시퀀스 가정(모든 스탠드가 같은 궤적을 따른다) 없이, 측정된 rate를 쓴다.
- "aggradation이 임령과 함께 둔화된다"는 정성적으로는 교과서(Bormann & Likens 1979; Franklin et al. 2002)이지만, **대륙 규모에서 직접 측정된 rate로, 지표별 타임스케일을 분리해, 기능(GPP)과 연결한** 정량화는 드물다.

## 핵심 질문 / 가설
1. **Q1 (rate–age scaling):** 구조 발달 속도(LAI·height·VCI·gap·heterogeneity 추세)는 임령에 따라 감쇠하는가? 함수 형태는(지수 감쇠? 임계 전이?).
   - H1: 초기 aggradation(height·LAI·VCI 상승) 속도는 임령↑에 따라 감소.
2. **Q2 (metric-specific timescales):** 지표마다 성숙 타임스케일이 다른가?
   - H2: 높이/LAI는 먼저 포화, **gap dynamics(deep-gap 추세)는 노령림에서 뒤늦게 증가**(aggradation→gap-phase 전이).
3. **Q3 (structure–function coupling):** 구조 aggradation 속도가 **GPP 추세**와 결합되는가?
   - H3: 구조 축적이 빠른 젊은 스탠드일수록 GPP 추세도 양(+); 노령림은 둘 다 정체.

## 데이터 (이미 확보)
- **반복 구조**: `data/per_year_v2_26.csv` — 26사이트 per-plot per-year 구조 13지표(Canopy_Ht, Max_Ht, Rumple, Rugosity, Deep_Gap, Vert_SD, Vert_CV, Gini, VCI, FHD, LAI, Q95, Ht_Ratio). 파생 추세: `lidar_pooled_predictors_26.csv`의 `*_trend`, `*_nyears`.
- **임령**: `data/plot_stand_age_gami_26.csv` (GAMI v3.1, 100m, 2020, 20-member ensemble mean; 원본 대비 r=0.83 검증).
- **기능**: `data/plot_pml_gpp_ts_26.csv` (PML-V2 2000–2024 연별 + GPP trend/sd).
- **교란**: `data/plot_disturbance_neon_26.csv` (NEON DP1.10111 + plantStatus) — 교란 스탠드 제외/통제용.

## 분석 계획
1. **추세 추정**: 플롯별 지표 연간 slope (관측 ≥3회), 노이즈 대비 효과크기 산출.
2. **rate–age**: slope ~ age (도메인 FE + 사이트 클러스터 SE; 또는 베이지안 multilevel). 함수형(선형/지수/조각별) 비교(AIC).
3. **지표별 타임스케일**: 각 지표의 rate–age 곡선에서 "성숙 나이"(rate가 0에 근접하는 임령) 추정·비교.
4. **구조–GPP 결합**: 구조 rate ↔ GPP trend, 임령 조건부.
5. **강건성(아래 리스크 §의 핵심)**: 획득 횟수·point density·acquisition 조건 민감도, 교란 스탠드 제외 재적합, GAMI 대체 임령(30m) 재현.

## 리스크 & 완화 (★ = killer)
- **★ 추세 신뢰성**: 반복 창이 ~5–8년, 플롯당 3–6회. "추세"가 생물학적 변화인지 **acquisition 간 불일치(BRDF·point density·registration·phenology) 아티팩트**인지 먼저 못 박아야 함. → 노이즈 플로어 대비 효과크기, n_years·point-density 민감도, 동일 사이트 내 재현성, null(무변화) 시뮬레이션.
- **임령 유효 n**: 임령이 사실상 사이트수준(ICC 높음) → 유효 n≈26, 검정력 제약. → 베이지안 partial pooling, 사이트수준 요약 병행.
- **"known result" 리젝**: 단순 "감쇠"면 기각. → Q2(지표별 타임스케일) + Q3(GPP 결합)로 정량·다차원·기능 연결이 반드시 필요.
- **임령·생물지리·수종 교란(confound)**: 임령이 도메인/수종과 상관. → 도메인 FE, 수종·교란 통제, partial regression.

## 타깃 저널
- 1순위: **Global Change Biology** / **Agricultural and Forest Meteorology** (구조–GPP–탄소 각도 살릴 때)
- 대안: **Forest Ecology and Management**, **Remote Sensing of Environment**(방법 각도)

## 다양성 논문(UNIFIED)과의 경계
- UNIFIED 논문에는 age→구조추세를 **한 줄 validity check**로만 남긴다(동역학→nestedness 신호가 노이즈 아님을 보증). 
- 본 논문은 **다양성을 다루지 않는다** — 순수 구조 발달·기능 축. 두 논문이 데이터는 공유하되 질문은 직교.

## Next steps
1. 추세 신뢰성 검증(★) 파일럿 — 노이즈 플로어 vs 효과크기.
2. Q1 rate–age 함수형 적합(13지표).
3. Q2 타임스케일 랭킹 그림 초안.
4. 파일럿 결과 보고 후 full/kill 결정.
