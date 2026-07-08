# NEON_v2 — 정리된 최종 분석 (단일 BRDF 출처)

> 2026-06-30 생성. 기존 분석(분광 VI 출처 혼재)을 발견·수정하고 **NEON .002 BRDF 단일 출처**로 재구성.
> 기존 모든 분석은 `../_ARCHIVE_v1/` 한 폴더에 보관됨.

## 무엇이 바뀌었나 (핵심 1줄)
분광 VI를 **연도별 혼재 출처**(일부 NEON-BRDF + 일부 비-BRDF 우리계산)에서
→ **NEON .002 bidirectional BRDF 단일 출처(8개 연도)** 로 통일. 결론은 동일. 상세: `docs/00_DATA_SOURCES.md`.

## 폴더 구조
```
NEON_v2/
├── README.md                  ← 이 파일
├── docs/
│   └── 00_DATA_SOURCES.md     ← ⭐ 변수별 데이터 출처 명세 (먼저 읽기)
├── data/
│   ├── plot_vi_neon_brdf.csv          신규: NEON .002 BRDF VI, 8연도, plot×year (1856행)
│   ├── lidar_pooled_predictors.csv    재사용: LiDAR mean/sd/trend + 다양성 응답 (BRDF 무관)
│   ├── plot_sampling_metrics.csv      재사용: Hill q1/q2, richness, coverage
│   ├── beta_coverage_filtered.csv     재사용: β/LCBD turnover·nestedness
│   ├── lcbd_rarefied.csv              재사용: rarefied LCBD
│   └── FINAL_v2_pooled.csv            ⭐ 최종 통합 (772 plots × 81; 분석 입력)
├── scripts/
│   ├── 01_extract_vi_neon_brdf.py     .002 zip→plot VI 직접 추출 (단일 출처 보장)
│   ├── 02_build_dataset.py            VI 시간특징 + LiDAR + 다양성 병합
│   └── 03_fit_models.py               혼합모델 + forest figure
├── results/
│   └── v2_coeff.csv                   197 모델 계수 (response×predictor×feature)
└── figures/
    └── Fig_v2_forest.png              Hill q1 ~ RS (mean/sd/trend), 구조 vs 분광
```

## 재현 (순서대로)
```bash
python NEON_v2/scripts/01_extract_vi_neon_brdf.py   # ~25분 (E:/ .002 zip 읽기)
python NEON_v2/scripts/02_build_dataset.py
python NEON_v2/scripts/03_fit_models.py
```

## 핵심 결과 (cov≥0.9, 525 plots, 19 sites)
- **알파(Hill q1/q2)**: 분광 SAVI/EVI **양(+0.38~0.41)**, 구조 Deep_Gap **음(−0.39)**
- **Nestedness**: Gini **음(−0.31)**
- 전부 기존 결론과 일치 → 출처 정리는 방법론적 엄밀성만 개선

## 데이터 원칙
- 원본 절대 덮어쓰지 않음 (모두 신규 경로 출력)
- 분광=NEON BRDF 단일출처 / LiDAR·조성=BRDF 무관 전연도
- PRI·fPAR 제외 (cross-year 아티팩트)
