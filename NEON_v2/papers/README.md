# 논문 3편 구성 (NEON RS × 종다양성)

> NEON 19개 산림 사이트(CONUS). 각 폴더에 그림·결과 사본 정리(원본은 상위 `figures/`·`results/`). 캐논 총괄: `../docs/03_OBJECTIVES_COMPLETE.md`.

| 논문 | 축 | 한 줄 |
|---|---|---|
| **P1 — 패턴** | 다양성 예측 + 조절 | RS(구조 다양성 + 분광 VI)로 종다양성 예측, 교란·경관·천이가 그 관계를 조절 |
| **P2 — 프로세스** | 구조 발달 | 다중시기 LiDAR로 본 산림 구조변화 — 천이(임령) + 교란(press/pulse) |
| **P3 — 생산성** | 생산성–다양성 | 위성 DHI로 생산성–다양성 혹형 검출 (biome 내부에도 견고) |

## 공통 데이터
- 통합: `../data/FINAL_v2_full.csv`, 분광 = NEON .002 BRDF 개별 VI(NDVI·EVI·ARVI·SAVI, PRI·fPAR 제외)
- 임령: `../data/plot_stand_age_{gami,30m}.csv` (GAMI 100m = 유효 임령, ICC=0.75 사이트단위)

## 논문 간 대칭 구조
교란·임령이 **조성(P1)과 구조(P2)에 각각** 작용:
- 교란 → 조성(희귀종 손실, P1) / 구조(press·pulse, P2)
- 임령 → RS–다양성 관계 조절(P1) / 구조 발달(P2)
- 생산성 → 다양성 혹형(P3, between-site)

## 방법론적 핵심 (reviewer 대응)
- **겉보기 혹형 검증 필수**: 임령–다양성 혹형 = **생물지리 아티팩트**(forest type 통제 시 소멸, p=0.98) / DHI–다양성 혹형 = **진짜**(forest type 통제에도 p<1e-5, biome 내부 존재). → 위성 혹형은 반드시 forest type 통제.
- 분광 = 개별 VI(합친 분광다양성 아님), 임령 효과는 사이트 수준(유효 N≈19).
