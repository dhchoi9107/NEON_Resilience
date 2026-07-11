# 천이·생장 축 — 시간에 따른 복잡화와 임령(stand age)

> 2026-07-09. 가설: **시간↑ → 식생 생장 → 구조 복잡도↑.** 이를 (A) 시간 trend의 정보량, (B) 임령(succession) 두 축으로 검증.
> **귀속**: (A) 시간 trend 증분 → **Obj1**(RS 피처가 다양성 예측에 기여하나, 다리). (B) 임령 조절·복잡화 차이 → **Obj2 축③ 천이**(관계의 시간 맥락 조절자, recency와 동일 메커니즘의 장기 버전).

## 배경 — 데이터가 생장 가정을 지지
미교란 plot에서 구조 복잡도가 방향성 있게 증가: **FHD·VCI 85%↑, Max_Ht 82%↑, Gini 77%↑, Deep_Gap 76%↓(갭 닫힘)**. 평균 수고·Rumple은 혼재(갭 메우며 저층 식생 증가 + 8년은 짧음). → "생장→복잡화"는 특히 **수직 복잡도**에서 실재.

---
## (A) 시간 trend는 mean 대비 예측력을 더하는가? — YES, 구조만
`div ~ z(mean) + z(trend) + C(domain) + (1|site)`, trend 계수의 Wald 검정 = 증분 기여. (`130_trend_increment.py`)

- **60쌍 중 19쌍에서 trend가 mean 통제 후에도 유의. 19개 전부 구조(Structural), 분광 VI는 0개.**
- 응답별: 알파 5, **Turnover(베타) 7, Nestedness(베타) 7** → trend 정보는 **베타 다양성(조성 전환)** 에 집중.
- **해석**: 분광 VI의 시간변화는 정보 없음(site-mean 엽록소가 전부, static). **구조의 생장 궤적(trend)은 mean이 못 담는 조성 정보를 추가** — 특히 베타. 생장→복잡화가 베타 다양성에서 실측 신호.
- 그림 `M01_trend_increment.png`, 결과 `trend_increment.csv`.

---
## (B) 임령(stand age) — 천이단계 축

### 데이터 (외부, plot lon/lat 샘플)
NEON은 plot 임령 미제공 → 외부 산림연령 지도 2종 교차검증:
| 출처 | 해상도 | 방식 | 추출 |
|---|---|---|---|
| **GAMI v2.1** (Besnard 2024) | 100m | 전지구 NetCDF(218GB), FIA 4만+ XGBoost, 2020 | ✅ h5py 바이트레인지 lazy 추출(전체 다운 없이), 676/797 plot, 57초 |
| **Global 30m** (Zhang 2025) | 30m | Landsat CCDC 변화점, 자연/조림 | ✅ 북미 타일 11.4GB(NF 542+PF 515 tif) → 709/797 plot |

### 교차검증 (`stand_age_crossval.csv`, 그림 `N02_age_crossval.png`)
- **Global 30m은 40년에서 포화**(76% plot이 ~40) — Landsat 1985–2024 레코드 길이가 상한. 1985 이후 무교란 성숙림은 전부 캡. **진짜 임령이 아니라 "Landsat 시대 내 마지막 교란 이후 경과년(≤40)".**
- 두 제품 일치: 전체 Spearman ρ=0.16(포화 탓 약함), **젊은층(GAMI<40, n=82) ρ=0.57(p=3e-8)**, 젊음(<30) 분류 82% 일치.
- WREF·TEAK·RMNP(GAMI 170–204년, 노령림)는 30m=40 포화. → **성숙림 위주 NEON엔 GAMI가 유효 임령**, 30m은 GAMI 젊은 추정치를 독립 검증 + "최근교란 recency" proxy(Obj2 연결).

**GAMI 사이트 중앙 임령 — 생태적으로 검증됨:**
WREF **202년**(old-growth 미송 ✔), TEAK 204·RMNP 170·SOAP 157(서부 노령 침엽수), ABBY **14년**(재생 관리림 ✔), 동부 성숙림 BART 108·HARV 101·SCBI 106, JERC 31·OSBS 32(FL 장엽송, 빈번 산불로 리셋 ✔).

### 결과 (GAMI, cov≥0.9 & 임령 있는 471 plot; `140_stand_age.py`)

**(i) 다양성 ~ 임령: 선형효과는 없으나 혹형(∩)이 존재 — 단 생물지리 교란** (`N09`).
- 선형은 null(Hill q1 lin p=0.54). 그러나 **2차 혹형(∩)**: Hill q1/q2·turnover가 **중간연령에서 최대** — plot p<1e-10, **사이트평균(n=19) p<0.01**. **nestedness만 반대 U자(∪)**(p=0.009; 젊음·노령서 높고 중간서 낮음 — 희귀종 nested 패턴과 정합).
- **forest type 통제 시 완전 소멸**(`N10`·`146_foresttype_control.py`, `foresttype_hump.csv`): 혹형 2차 p가 통제無 4.8e-13 → **+forest_type(문서화 3범주) 0.98 / +상록성(DHI_min 연속) 0.95**. **domain·forest type·상록성 어느 것으로 통제해도 혹형 소멸** → **생물지리(forest type) 교란 확정**. 중간연령 ≈ 동부 낙엽수림(고 수종다양성 ~85–110년), 양끝 ≈ 저다양성(젊은 SE 소나무·관리림 + 노령 서부 침엽수). **type 내부엔 혹형 없음**(낙엽림은 나이들수록 오히려 감소, 상록림 평평). **순수 임령혹형 아님.**
- ↔ 대조: **Obj3 DHI 생산성혹형은 통제 후에도 견고**(다항식 p<0.001, GAM EDF 3.2). 임령혹형은 교란 소멸이라 약함.
- 검증: DHI_min(상록성 대리) evergreen 0.36 > deciduous 0.29 > mixed 0.26로 방향 일치.
> ⚠ **임령의 스케일 한계**: 임령 분산의 **75%가 사이트 간(ICC=0.75)**, within-site 25%는 대부분 GAMI 100m 아티팩트(ABBY 중앙 14년인데 내부 SD 70년; 성숙림은 내부 SD ~7년으로 거의 단일연령). → 임령은 **사이트 단위 변수**(유효 N≈19 사이트/9 도메인). **site 통제 후 남는 임령잔차(`N08`)는 노이즈라 의미 없음**(데이터 한계). 임령 효과·조절은 plot-level 유의성 아닌 **사이트 수준**으로 보수적 해석(의사반복 주의).
> ⚠ **GAMI 정확도**: v2.1(100m)은 **독립 검증 미보고**(Besnard 2021 상속: NSE 0.60, **RMSE ~48년**, 젊음 과대·노령 과소추정 편향). 제작자 권고 = **최소 2-decade(20년) binning**. → **20년 구간 재실행(`N11`, `age_binned_*.csv`)에서 결과 견고**: 구조변화 4지표(VCI −0.002 p=5e-4, LAI, Deep_Gap +0.004 p=5e-5, Canopy_Ht) 연속과 동일, 조절효과 10/11 유의 전부 음. → **절대 임령 아닌 20년 구간으로 해석하면 편향에도 결론 유지.**

**(ii) 임령이 RS–다양성 관계를 조절 — 11/19 유의, 전부 음(−):**
| 관계 | 상호작용 | p |
|---|---|---|
| Hill q1 ~ LAI × age | −0.214 | 0.002 |
| Hill q1 ~ VCI × age | −0.168 | 0.003 |
| Hill q1 ~ SAVI × age | −0.179 | 0.005 |
| Hill q2 ~ EVI × age | −0.178 | 0.007 |
| Turnover ~ EVI × age | −0.146 | 0.012 |
- **핵심**: RS(분광 SAVI/EVI/NDVI + 구조 LAI/VCI/FHD)–다양성 커플링은 **젊은 임분에서 강하고 성숙할수록 약해짐**. 성숙 폐쇄 캐노피에서 분광·구조 신호가 포화→종다양성과 디커플링. **님 생장 가설의 핵심 귀결.**
- ⚠ **단, 임령은 사이트 단위(ICC=0.75)** → 이 조절은 본질적으로 **사이트 수준 현상**(유효 N≈19 사이트/9 도메인). plot-level p는 낙관적일 수 있어 **사이트 수준 조절**로 보수적 해석(의사반복 주의). 방향·효과크기는 일관되나 강한 단정은 회피.

**(iii) 복잡화 속도(구조 trend) ~ 임령 — 4/8 유의, 천이 교과서:**
| trend | β(age) | p | 의미 |
|---|---|---|---|
| VCI_trend | **−0.248** | 0.0009 | 젊을수록 수직복잡도 빨리 축적 |
| LAI_trend | **−0.218** | 0.003 | 젊을수록 엽면적 빨리 축적 |
| Deep_Gap_trend | +0.283 | 0.0002 | 노령림은 갭 열림(gap-phase 동태) |
| Canopy_Ht_trend | +0.259 | 0.001 | 노령림 수고 trend 양 |
- **핵심**: 젊은 임분은 엽면적·수직복잡도를 빠르게 쌓고(생장), 노령림은 갭-단계 동태(교란·고사로 갭 열림). 생장→복잡화 가정을 지지하며 **복잡화의 *종류*가 천이단계에 따라 전환**됨을 보임.
- **산점도 `N03`** (점 색=site, bivariate): VCI_trend **r=−0.41**, FHD_trend **r=−0.40**(젊을수록 수직복잡도·엽층다양성 빠르게 축적→~200년 정체), Deep_Gap_trend **r=+0.29**(노령 갭 열림), LAI/Rugosity/Max_Ht/Canopy_Ht 음, Gini ns.
  - **site 색이 드러내는 것**: 임령 구배가 부분적으로 **사이트 간** 구배(ABBY ~15년 복잡화 최고 ↔ WREF/TEAK ~200년 정체 ↔ 동부 HARV/MLBS/SCBI 중간 ~90–110년). → 순수 임령효과는 site 랜덤효과 통제한 혼합모델(위 표·N01)이 더 보수적·타당.
- **부분회귀 `N04`** (도메인·site 랜덤효과 통제 후 partial residual = 순수 임령효과, 140 모델 계수와 일치):
  - **살아남는 4개**: Deep_Gap_trend **+0.28***, Canopy_Ht_trend **+0.26**, VCI_trend **−0.25***, LAI_trend **−0.22** → 천이단계별 복잡화 전환이 site 넘어 실재(젊음=수직복잡도·엽면적 축적, 노령=갭 열림+생존목 수고 성장).
  - **소멸**: **FHD_trend는 bivariate r=−0.40 → 통제 후 −0.05 ns** = 사이트/도메인 교란(동·서부 차이)이었음. Rugosity·Max_Ht·Gini도 ns.
  - Canopy_Ht_trend는 bivariate 음→통제 후 양 **부호 반전**(Simpson's paradox). 분해(`N06`): **between-site r=−0.21**(젊은 ABBY 평균수고 급성장 +0.47 ↔ 노령 SOAP/TEAK 정체·감소) vs **within-site r=+0.16**(p=0.001, 사이트 내 성숙 patch 우점목 성장). 평균수고는 하층 재생에 민감(Max_Ht와 달리) → 사이트 간엔 "어린 숲 급성장", 사이트 내엔 "성숙 patch 성장"이 지배해 방향이 갈림. 부분회귀가 between 제거→within(양) 노출.
- **비교 패널 `N05`** (위=bivariate, 아래=부분회귀, 점·묶음=9개 도메인 convex hull): 도메인 hull이 소멸/견고를 시각화 — FHD는 젊은 도메인(D03 좌상)↔노령 도메인(D16·D17 우하) **도메인 간 배치**가 bivariate 기울기를 만들고 통제 시 평평(β=−0.05 ns). VCI(−0.25)·Deep_Gap(+0.28)은 도메인 *내*에서도 기울기 유지=견고. Canopy_Ht 반전.

### 30m 재실행 (robustness, 검열 감안)
30m(≤40 검열)으로 (B) 재실행 시 조절효과 패턴이 GAMI와 갈림(일부 반대 부호). **모순 아님** — 30m은 76% 포화라 사실상 "최근교란된 젊은 25%" 대비이므로 진짜 임령이 아닌 recency 축을 측정. 따라서 **(B)의 유효 결론은 GAMI(진짜 임령) 기준**이며, 30m은 젊은층 검증·recency proxy로만 사용. (Obj2 severity/recency 결과와 정합)

### 다양성 오버레이 (`N07` raw / `N08` partial)
- **`N07`(raw, bivariate)**: 임령 × VCI_trend, 점 크기 = 종다양성 4지수(색=site). 알파(Hill) 큰 점=중간연령·고 VCI_trend 동부 군집(MLBS/HARV), **nestedness 큰 점=젊은 ABBY**(재생 임분=nested 종빈약, 생태적 타당). ※ raw라 site 교란 포함.
- **`N08`(partial, added-variable)**: x·y·다양성 **모두 site 고정효과로 잔차화**. 크기=|다양성 잔차|, 색=부호(빨강 기대↑/파랑 기대↓). **부분 임령→VCI_trend r=−0.16***(site 통제 후에도 생존=견고). 다양성 잔차(빨강/파랑)는 공간에 **비교적 섞임** → site 걷어내면 임령×복잡화가 다양성 편차를 크게 조직하지 않음 = **임령 직접효과 없음(조절만)** 과 정합. 다양성 신호 대부분은 between-site(사이트 정체성).

### 종합 (한 문장)
> 임령은 다양성 *수준*을 직접 정하진 않지만, **RS–종다양성 관계의 강도(젊을수록 강)와 구조 복잡화의 양상(젊음=축적, 노령=갭동태)을 좌우**한다. 즉 원격탐사로 다양성을 추론할 때 **천이단계가 조절자**다.

그림 `N01_stand_age.png`. 결과 `stand_age_{models,moderation,trendage}_gami.csv`.

## 재현
```
python scripts/obj23/130_trend_increment.py          # (A)
python scripts/obj23/141_extract_age_30m.py          # 30m 타일 추출 (다운로드 후)
python scripts/obj23/140_stand_age.py                # (B) 혼합모델, 존재하는 age 소스 자동 병합
python scripts/obj23/142_stand_age_figures.py        # N02~N06 교차검증·산점도·부분회귀·Simpson
python scripts/obj23/143_age_diversity_bubble.py     # N07 임령×구조변화, 점크기=다양성 4지수
python scripts/obj23/144_age_diversity_partial_bubble.py # N08 partial 잔차화 버블
python scripts/obj23/145_age_diversity_scatter.py    # N09 임령 vs 4다양성지수 산점도
```
