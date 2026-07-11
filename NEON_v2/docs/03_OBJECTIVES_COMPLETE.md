# NASA NPP 제안서 — 3개 목적 총괄 (canonical)

> NEON 19개 산림 사이트(CONUS), 2013–2026. 통합 데이터 `data/FINAL_v2_full.csv`.
> 분광 = NEON .002 bidirectional BRDF 단일출처 개별 VI(NDVI·EVI·ARVI·SAVI), PRI·fPAR 제외.
> 최종 재정리 2026-07-10. 상세 목적별 문서: `01_OBJ23.md`(방법), `02_OBJ2_SYNTHESIS.md`(Obj2), `04_SUCCESSION_STANDAGE.md`(천이).

**한 눈에 — 세 목적의 논리:**
- **Obj 1 = 관계**: RS 다양성(구조+분광)이 종다양성을 예측하는가? (정적, 구조 vs 분광, 알파 vs 베타, 다중스케일)
- **Obj 2 = 조절**: 교란·경관 이질성·천이(임령)가 그 관계와 다양성을 시간에 따라 어떻게 바꾸는가?
- **Obj 3 = 생산성**: RS 다양성 ↔ 생산성(Sentinel DHI) — 생산성–다양성 혹형.

---
## Obj 1 — RS 다양성 ↔ 종다양성 (관계) ✅
**RS 예측변수(구조 LiDAR + 분광 개별 VI) ↔ 종다양성(알파 Hill q1/q2 + 베타 LCBD turnover·nestedness)**, 다중스케일 domain>site>plot.

| 응답 | 주요 예측변수 (mean, p<0.05) |
|---|---|
| Hill q1/q2 (알파) | 분광: **SAVI +0.38, EVI +0.38**, NDVI +0.28, ARVI +0.27 / 구조: LAI +0.36, VCI +0.31, Deep_Gap −0.39, Rumple −0.30, Gini +0.17 |
| LCBD turnover (베타) | 구조: VCI +0.14, Gini +0.12, Deep_Gap −0.15, Rumple −0.14 / **분광 VI 전부 비유의** |
| LCBD nestedness (베타) | 구조: Gini −0.31, Vert_CV −0.30, VCI −0.26, Canopy_Ht +0.24 / **분광 VI 전부 비유의** |

**핵심:**
- **알파 = 분광 VI(SAVI/EVI)가 최상위**(구조 LAI/VCI 동급·상회) + 구조. **베타 = 구조 전용**(분광 비유의).
- → 분광 VI = 알파(엽록소·생산성 신호), 구조 = 알파+베타(공간 조성 전환). 상호 보완.
- 다중스케일: 분산 대부분 domain>site>plot (F04).
- **다리→Obj2 (A) 시간 trend 증분**: RS 시간변화(trend)가 mean 대비 예측력 더하는 19/60쌍 **전부 구조(분광 0), 베타 집중** → 구조 생장궤적이 조성 정보 추가. 이 시간신호의 *맥락 조절*은 Obj2에서.

**결과 맵핑** — 스크립트 `01·02·03·04·05`, `obj23/110`(K01), `130`(A) / 그림 `F01`(forest)·`F02`(scatter)·`F03`(feature)·`F04`(분산분해)·`F05·F06`(도메인·사이트)·`F07·F08`(시계열·trend분포)·`F09`(커플링)·`Fig_v2_forest`·**K01**·**M01**(A) / 결과 `v2_coeff.csv`·`v2_variance_decomp.csv`·`trend_increment.csv`
> 폐기: 합친 분광다양성(Rao Q/FEve, `100_spectral_diversity.py`·`obj1_specdiv.csv`)은 Obj1에서 제외·보관.

---
## Obj 2 — 교란·경관·천이가 관계를 조절 (시간·맥락) ✅
관계(Obj1)를 **세 맥락축**이 시간에 따라 조절. 상세 `02_OBJ2_SYNTHESIS.md`.

### 축 ① 교란 (사건)
- **NEON 자체 기록 필수**: 원격(MTBS·Hansen·NBR)은 곤충 0/10·선택벌채 누락 → NEON 이벤트(DP1.10111)+현장 plantStatus로 정의.
- **BACI 시간분해**: 교란 → **희귀종 ~1종 손실**(8/8 사양, MannW p≤0.001), **우점 다양성(Hill)은 회복력**. 구조 회복력 ≠ 조성 회복력.
- **press/pulse**: press(곤충)=ΔLAI −0.18(고사), pulse(산불/풍해)=ΔLAI +0.22(재생). 다중시기 LiDAR가 체제 구분.
- **RTM 통제**: naive 복잡도완충은 대부분 평균회귀. 통제 후 **Deep_Gap(완충)·Gini(증폭)만 진짜**, 수고·LAI·FHD·richness는 RTM.
- **severity/recency 조절**: 분광–알파 커플링은 **갓 교란서 강→회복(recency↑)하며 약화**.
- 맵핑 — 스크립트 `30·31·40·62`(교란원)·`80·81·82`(NEON native)·`50·51`(모델)·`52`(severity/recency)·`70·71·72`(BACI)·`90·91·92`(Choi2023/RTM) / 그림 `G03·G04·G05`·`H01·H02·H03`·`I01·I02·I03`·`J01~J06` / 결과 `obj23_disturbance`·`obj2_severity_recency`·`baci`·`baci_sensitivity`·`choi2023`·`complexity_rtm`

### 축 ② 토지이용 이질성 (경관 맥락)
- **직접효과**: 파편화(edge density) → 다양성 **−0.30~−0.35**, 산림비율 +0.22~0.28.
- **조절효과**: **SAVI × 이질성 → turnover +0.21**(p<0.001) — 경관이 RS–종 관계를 바꿈.
- 맵핑 — 스크립트 `101` / 그림 **K02·K03** / 결과 `obj2_heterogeneity`

### 축 ③ 천이·임분연령 (시간 맥락) — 신규, ★임령은 여기
상세 `04_SUCCESSION_STANDAGE.md`. 임령 = **GAMI 100m**(진짜 임령 1–227, 218GB NetCDF lazy 추출) + **Global 30m**(≤40 포화=recency proxy, 교차검증 젊은층 ρ=0.57).
- **(i) 임령 직접효과: 선형 없음, 혹형(∩)은 있으나 생물지리 교란** — Hill/turnover 중간연령 최대(사이트평균 2차 p<0.01, nestedness는 U자∪), 단 domain 통제 시 소멸(중간연령=동부 고다양성 낙엽수). 순수 임령효과 아님(Obj3 DHI혹형은 통제 후 견고와 대조).
- **(ii) 임령이 RS–다양성 커플링 조절 — 11/19 유의 전부 음**(SAVI·EVI·LAI·VCI × age): **젊은 임분서 커플링 강, 성숙할수록 캐노피 포화·디커플링.**
- **(iii) 복잡화 속도가 임령에 따라 다름**: 젊을수록 VCI·LAI trend 큼(빠른 축적), 노령림은 Deep_Gap trend 양(갭-단계 동태).
- **Obj2 통합 논리**: 임령 조절 = 축①의 **recency×RS 조절과 동일 메커니즘의 다른 시간척도**(교란 recency 수년 ↔ 천이 임령 수십~수백년). 30m 임령이 recency로 포화되는 것이 이를 방증.
- 맵핑 — 스크립트 `140·141`(+`130` A는 Obj1) / 그림 **N01·N02** / 결과 `stand_age_{models,moderation,trendage}_{gami,30m}`·`stand_age_crossval`

**Obj2 종합**: 관계(Obj1)는 고정이 아니라 **교란(사건)·경관 파편화(공간)·천이단계(시간)** 세 맥락에 의해 조절된다.

---
## Obj 3 — RS 다양성 ↔ 생산성 (DHI) ✅
Sentinel-2 **DHI**(누적/최소/계절변동, GEE 서버사이드) ~ 종다양성.
- ★ **생산성–다양성 혹형(hump)**: Hill q1·q2 ~ DHI누적 = **inverted-U**(다항식 p<0.001, GAM EDF 3.2 p<0.001). 중간 생산성에서 다양성 최대.
- **선형 +0.12는 착시**(혹형 상승부만). 평균·선형 모델은 혹형을 놓침 → 다항식/GAM 필수.
- ★ **혹형 견고성 검증**(`L04`, `dhi_hump_confound.csv`): **forest type 통제해도 생존**(Hill q1 p=9.4e-6, turnover 4.8e-14), 상록성 통제도 생존, **evergreen·deciduous biome 내부에도 혹형 존재** → **생물지리 아티팩트 아님**. ↔ **임령–다양성 혹형은 forest type 통제 시 즉사(p=0.98)**로 근본적으로 다름. 단 between-site 거시생태 스케일(site 랜덤 통제 시 소멸=이 관계의 고유 스케일).
- 교란강도·구조·분광다양성은 선형/단조(중간교란가설 미지지).
- 맵핑 — 스크립트 `10·20·60·61`(DHI/Sentinel)·`120·121·122`(비선형/GAM/혹형)·`147`(혹형 검증)·`50`(obj23모델) / 그림 `G01·G02`·`L01·L02`·**L03**(혹형 headline)·**L04**(견고성) / 결과 `obj23_dhi_coeff`·`nonlinear`·`gam`·`dhi_hump_confound`

---
## 전체 그림 인덱스
| 코드 | 내용 | 목적 |
|---|---|---|
| F01–F09 | 종다양성~RS(forest/scatter/분산분해/도메인·사이트/시계열/커플링) | Obj1 |
| **K01** | 종다양성~RS(분광 개별 VI + 구조) | Obj1 |
| **M01** | 시간 trend 증분 예측력(구조만) | Obj1(A) |
| G03·G04·G05 | 교란 지도·커플링·회복 | Obj2① |
| H01–H03 | severity·recency·조절 | Obj2① |
| I01–I03 | BACI 전후·dose·변화결합 | Obj2① |
| J01–J06 | press/pulse·복잡도완충·RTM·계산법 | Obj2① |
| **K02·K03** | 토지이용 이질성 직접·조절 | Obj2② |
| **N01–N06** | 임령~다양성·복잡화 / 교차검증 / 산점도 / 부분회귀 / 도메인묶음 비교 / Canopy_Ht Simpson 분해 | Obj2③ |
| **N07·N08** | 임령×구조변화 점크기=다양성 (raw / partial 잔차화) | Obj2③ |
| **N09** | 임령 vs 종다양성 4지수 산점도(plot+사이트평균) | Obj2③ |
| **N10** | forest type 통제=혹형 소멸(생물지리 교란 확정) | Obj2③ |
| G01·G02 | 다양성~DHI(선형) | Obj3 |
| L01·L02 | 비선형(다항식·GAM) | Obj3 |
| **L03** | ★ 생산성–다양성 혹형 | Obj3 |

## 핵심 데이터·스크립트
- 통합: `data/FINAL_v2_full.csv` / 임령: `plot_stand_age_{gami,30m}.csv`
- 분광 단일 BRDF 출처: `01·02` (`docs/00_DATA_SOURCES.md`)
- 스크립트: `scripts/`(Obj1 핵심 01~05), `scripts/obj23/`(Obj2·3 = 10~101, 천이 130·140·141, 그림 110~122)
