# P1 — 패턴: RS로 종다양성 예측 + 무엇이 그 관계를 조절하나

## 주제
원격탐사 유래 다양성(**구조 다양성 LiDAR + 분광 개별 VI**)으로 종다양성(알파 Hill q1/q2 + 베타 LCBD turnover·nestedness)을 다중스케일(domain>site>plot) 예측하고, **교란·경관 이질성·천이(임령)가 그 RS–종 관계를 어떻게 조절**하는지 규명.

## 핵심 결과
1. **알파 = 분광 VI(SAVI/EVI +0.38)가 최상위** 예측변수(구조 LAI/VCI 동급·상회), **베타 = 구조 전용**(분광 VI 비유의). 구조와 분광이 상호 보완.
2. **시간 trend 증분**: RS 시간변화가 mean 대비 예측력 더하는 19/60쌍 전부 구조(분광 0), 베타 집중.
3. **교란이 관계·조성 조절**: BACI로 희귀종 ~1종 손실(우점 다양성은 회복력); 분광–알파 커플링은 갓 교란서 강→recency↑ 약화.
4. **경관 이질성**: 파편화 → 다양성 −0.30~−0.35; SAVI × 이질성 → turnover +0.21(관계 조절).
5. **임령 조절**: RS–다양성 커플링 젊을수록 강(사이트 수준, 유효 N≈19). **단 임령 主효과는 null.**
6. ★ **방법론**: 임령–다양성 **혹형(∩)은 생물지리 아티팩트** — forest type/상록성 통제 시 소멸(p 4.8e-13→0.98). 위성 기반 혹형 해석 시 forest type 통제 필수.

## 그림 (figures/)
- **K01** 종다양성~RS(분광 개별 VI + 구조) ★headline
- F01–F09 forest/scatter/feature/분산분해/도메인·사이트/시계열/trend분포/커플링
- Fig_v2_forest (Hill q1 forest), **M01** 시간 trend 증분
- G03·G04·G05 교란 지도·커플링·회복 / H01·H02·H03 severity·recency·조절 / I01·I02·I03 BACI
- **K02·K03** 토지이용 이질성 직접·조절
- **N09** 임령~다양성 4지수 / **N10** forest type 통제=혹형 소멸(생물지리 교란)

## 결과 (results/)
v2_coeff, v2_variance_decomp, obj1_specdiv(보관), trend_increment, obj2_heterogeneity, obj2_severity_recency, baci, baci_sensitivity, obj23_disturbance, stand_age_moderation_{gami,30m}, stand_age_models_{gami,30m}, foresttype_hump, stand_age_crossval

## 관련 스크립트 (../../scripts/)
`03·04·05`, `obj23/110`(K01)·`130`(trend)·`50·51·52·70~72·80~82`(교란)·`101`(경관)·`145·146`(임령 다양성·혹형)
