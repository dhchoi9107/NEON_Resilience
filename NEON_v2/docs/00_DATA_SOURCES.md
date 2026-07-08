# NEON_v2 — 데이터 출처 명세 (단일 출처 원칙)

> 작성 2026-06-30. 기존 분석에서 **분광 VI가 연도별로 출처가 섞여** 있던 문제를 발견·수정한 기록.
> 핵심 원칙: **모든 분광 VI는 NEON .002 bidirectional(BRDF 보정) 단일 출처에서만 사용.**

---

## 1. 왜 다시 했나 — 발견한 문제

기존 `plot_spectral_1m`(= 분석에 쓰인 분광 데이터)은 **14개 연도(2013–2026)** 를 담고 있었으나,
NEON의 BRDF 보정 VI 제품(**DP3.30026.002 "bidirectional"**)은 **8개 연도만** 존재한다:

| | 연도 |
|---|---|
| NEON .002 BRDF VI **있음** | **2016, 2017, 2018, 2022, 2023, 2024, 2025, 2026** |
| NEON .002 BRDF VI **없음** | 2013, 2014, 2015, 2019, 2020, 2021 |

→ 빠진 6개 연도의 기존 값은 **우리가 raw 초분광에서 계산한 것**이었고, 우리 BRDF 보정은
비율 VI(NDVI 등)에 사실상 무효(no-op)였음 → **그 6개 연도는 사실상 비-BRDF**.
즉 분석 VI가 **연도별로 BRDF 처리가 일관되지 않았다.** (디스크·NEON API 양쪽 확인)

검증: 2022는 새 .002 추출과 기존이 **NDVI r=0.998, Δ=0.000**(동일 출처) /
2017은 **r=0.985, Δ=−0.032**(출처 차이 흔적).

---

## 2. 결정 (옵션 A 채택)

**NEON .002 BRDF VI가 존재하는 8개 연도만 사용.** 비-BRDF 6개 연도는 분광 분석에서 제외.
→ 분광 시간특징(mean/sd/trend)이 **전부 동일한 NEON-BRDF 기하**에서 계산됨.

---

## 3. 변수별 출처 (최종)

| 변수군 | 출처 | 연도 | BRDF | 비고 |
|---|---|---|---|---|
| **분광 VI** (NDVI, EVI, ARVI, SAVI) | **NEON DP3.30026.002 bidirectional** | 2016,17,18,22–26 (8) | ✅ NEON FlexBRDF+지형 | `plot_vi_neon_brdf.csv`, 신규 직접 추출 |
| ~~PRI, fPAR~~ | (제외) | — | — | cross-year 알고리즘 아티팩트 ([[07_spectral_calibration]]) |
| **LiDAR 구조** (FHD, LAI, VCI, Rumple, Deep_Gap, Gini, Canopy_Ht 등 13종) | NEON DP1.30003.001 LiDAR | 전 연도 | 해당없음(기하) | `lidar_pooled_predictors.csv`, 10m 직접 계산 |
| **분류 다양성 응답** (richness, Hill q1/q2) | NEON 식생구조 DP1.10098.001 | pooled(전 연도) | 해당없음 | `plot_sampling_metrics.csv`, coverage≥0.9 |
| **β/LCBD 응답** (turnover/nestedness rare) | 위 조성에서 Baselga 분해 | pooled | 해당없음 | `beta_coverage_filtered.csv`, `lcbd_rarefied.csv` |

핵심: **BRDF가 의미 있는 건 분광뿐.** LiDAR(기하)·분류조성(현장조사)은 BRDF와 무관 → 그대로 재사용.

---

## 4. 분광 VI 연도 분포 (8 BRDF 연도, 661 plots)

| nyears | plots |
|---|---|
| 4년 | 129 |
| 3년 | 279 |
| 2년 | 250 |
| 1년 | 3 |

mean = 가용 연도 평균, sd = 연간 변동(≥2년), trend = OLS 기울기(≥3년).

---

## 5. 결과: 출처 정리해도 결론 불변

`response ~ z(predictor) + C(domain) + (1|site)`, coverage≥0.9, 525 plots / 19 sites / 197 models.

| 응답 | 최상위 예측변수 | 기존 대비 |
|---|---|---|
| Hill q1 (mean) | **SAVI +0.38, EVI +0.38**, NDVI +0.28, ARVI +0.27 / **Deep_Gap −0.39** | 기존 SAVI+0.39·EVI+0.38·Deep_Gap−0.37 **일치** |
| Hill q2 (mean) | SAVI +0.41, EVI +0.40 | 일치 |
| Nestedness (mean) | Gini −0.31 | 일치 |

→ **단일 NEON-BRDF 출처로 재계산해도 과학적 결론 동일.** 방법론만 깨끗해짐.

## 6. 논문 기술 문구(안)
> "Spectral vegetation indices were taken solely from NEON's BRDF- and topographically-corrected
> bidirectional product (DP3.30026.002), available for 8 acquisition years (2016–2018, 2022–2026).
> Years lacking the bidirectional product were excluded to ensure a single, consistent BRDF-corrected
> source. PRI and fPAR were excluded due to a cross-year calibration artifact. LiDAR structural metrics
> and field-based taxonomic composition are BRDF-independent and use all available years."

관련: [[07_spectral_calibration]] (fPAR/PRI 제외), BRDF near-nadir 검증(모자이크 관측각 중앙 ~5°)
