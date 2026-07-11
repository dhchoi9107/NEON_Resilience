# P3 — 생산성: 위성 DHI로 본 생산성–다양성 혹형

## 주제
Sentinel-2 **DHI(Dynamic Habitat Indices: 누적/최소/계절변동, GEE 서버사이드)** ~ 종다양성. 위성 생산성 대리변수로 **생산성–다양성 혹형(hump)**을 CONUS 산림 전역에서 검출·검증.

## 핵심 결과
1. ★ **생산성–다양성 혹형(inverted-U)**: Hill q1·q2 ~ DHI누적 = **∩** (다항식 p<0.001, GAM EDF 3.2 p<0.001). 중간 생산성에서 다양성 최대. **선형 +0.12는 착시**(혹형 상승부만).
2. **혹형은 진짜 — biome 아티팩트 아님** (`L04`, `dhi_hump_confound.csv`): forest type 통제해도 **생존(Hill q1 p=9.4e-6)**, 상록성 통제도 생존(2.4e-5), **evergreen·deciduous biome 내부에도 혹형 존재.**
   - ↔ **대조**: 임령–다양성 혹형은 forest type 통제 시 즉사(p=0.98) = 순수 생물지리 아티팩트. **DHI 혹형은 이와 근본적으로 다름.**
3. **스케일**: between-site 거시생태 관계(site 랜덤 통제 시 p=0.12로 사라짐 — 이는 관계 자체가 사이트 간 스케일이기 때문, PDR의 고전적 스케일). 한계: 이 스케일에선 생산성이 기후·biome과 얽힘(모든 PDR 연구 공통).
4. minimum/variation DHI는 비유의; 교란강도·구조는 선형/단조(중간교란가설 미지지).

## 그림 (figures/)
- **L03** 생산성–다양성 혹형 ★headline / L01·L02 비선형(다항식·GAM)
- **L04** 혹형 견고성 검증(forest type별 곡선 — biome 내부에도 혹형)
- G01 DHI forest / G02 DHI 산점도

## 결과 (results/)
obj23_dhi_coeff, nonlinear, gam, **dhi_hump_confound**(통제별 2차 p)

## 데이터
`../../data/plot_dhi_sentinel.csv` (dhi_cum/min/var, Sentinel-2 2016+)

## 관련 스크립트 (../../scripts/obj23/)
`10·20·60·61`(DHI/Sentinel GEE)·`120·121·122`(비선형·GAM·혹형)·`147`(혹형 교란 검증 L04)
