# NASA NPP 제안서 — 3개 목적 완료 총괄

> NEON 19개 산림 사이트. 통합 데이터 `FINAL_v2_full.csv` (종다양성+구조다양성+분광지수+분광다양성+DHI+교란+토지이질성).
> 최종 점검 2026-07-01. 세 목적 모두 데이터 검증 완료.

---
## Obj 1 — RS 다양성 ↔ 종다양성 (다중스케일) ✅
**RS 예측변수(구조 LiDAR + 분광 개별 VI) ↔ 종다양성(알파 Hill + 베타 LCBD)**, domain>site>plot.
> 분광 = 개별 VI(NDVI·EVI·ARVI·SAVI) 직접 사용. 합친 분광다양성(Rao Q·FEve·FDiv)은 Obj1에서 제외(script 100·obj1_specdiv.csv는 보관용). PRI·fPAR 제외 유지.

| 응답 | 주요 예측변수 (mean, p<0.05) |
|---|---|
| Hill q1/q2 (알파) | 분광: **SAVI +0.38, EVI +0.38**, NDVI +0.28, ARVI +0.27 / 구조: LAI +0.36, VCI +0.31, Deep_Gap −0.39, Rumple −0.30, Gini +0.17 |
| LCBD turnover (베타) | 구조: VCI +0.14, Gini +0.12, Deep_Gap −0.15, Rumple −0.14 / **분광 VI 전부 비유의** |
| LCBD nestedness (베타) | 구조: Gini −0.31, Vert_CV −0.30, VCI −0.26, Canopy_Ht +0.24 / **분광 VI 전부 비유의** |

- **알파: 분광 개별 VI(SAVI/EVI)가 최상위 예측변수**(구조 LAI/VCI와 동급·상회), **베타: 구조 전용**(분광 VI 비유의)
- → 분광 VI = 알파(엽록소·생산성 신호), 구조 = 알파+베타(공간 조성 전환) 담당. 상호 보완.
- 그림 **K01** / 스크립트 03·110 / 결과 v2_coeff.csv

## Obj 2 — 교란 + 토지이용 이질성 (시간) ✅
상세: `02_OBJ2_SYNTHESIS.md`. 핵심:
- **교란**: NEON 자체 기록(이벤트+plantStatus) 필수(원격은 곤충 0/10 누락). BACI로 시간분해 →
  **희귀종 ~1종 손실, 우점구조 회복력**. Choi2023 press/pulse. RTM 통제 후 **Deep_Gap·Gini만 진짜 완충**.
- **토지이용 이질성**: **파편화 → 다양성 −0.35**, 산림비율 +0.28. **이질성이 RS–종 관계 조절**(SAVI×het +0.21).
- **핵심**: 구조적 회복력 ≠ 조성적 회복력 + 경관 파편화가 다양성·관계를 좌우.
- 그림 G03·I01·I02·I03·J01·J05·J06·**K02·K03** / 스크립트 50·70~72·80~82·90~92·101

## Obj 3 — RS 다양성 ↔ 생산성 (DHI) ✅
- Sentinel-2 **DHI**(누적/최소/계절변동, GEE 서버사이드) ~ 종다양성.
- ★ **생산성–다양성 혹형(hump)**: Hill q1·q2 ~ DHI누적 = **inverted-U** (다항식 p<0.001, GAM EDF 3.2 p<0.001).
  중간 생산성에서 다양성 최대 — 고전 생산성-다양성 관계를 Sentinel DHI로 검출. (선형 +0.12는 착시)
- **비선형이 핵심**: 평균·선형 모델은 혹형을 놓침. 교란강도·구조는 선형(중간교란가설 미지지).
- 그림 **L03**(혹형 headline)·L01·L02·G01·G02 / 스크립트 60·61·120·121·122 / 결과 gam.csv·nonlinear.csv

---
## 전체 그림 인덱스
| 코드 | 내용 | 목적 |
|---|---|---|
| F01–F09 | 종다양성~RS(forest/scatter/분산분해/시계열/커플링) | Obj1 |
| **K01** | 종다양성~RS예측변수(분광 개별 VI + 구조) | Obj1 |
| G01–G02 | 다양성~Sentinel DHI (선형) | Obj3 |
| **L03** | ★ 생산성–다양성 **혹형**(GAM) | Obj3 |
| L01–L02 | 비선형(다항식·GAM) 전체 검정 | Obj1–3 |
| G03 | 교란 지도 | Obj2 |
| I01–I03 | BACI(전후 다양성변화·회복·변화결합) | Obj2 |
| J01 | press/pulse 구조변화 | Obj2 |
| J05–J06 | 복잡도 완충·RTM 통제·계산법 | Obj2 |
| **K02–K03** | 토지이용 이질성 직접·조절효과 | Obj2 |

## 핵심 데이터·스크립트
- 통합: `data/FINAL_v2_full.csv`
- 분광 단일 BRDF 출처: 01·02 (docs/00_DATA_SOURCES.md)
- 스크립트: `scripts/` (Obj1 핵심), `scripts/obj23/` (Obj2·3, 50~110)
