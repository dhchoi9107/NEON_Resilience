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

**(i) 다양성 ~ 임령: 직접효과 없음** (Hill q1 선형 p=0.54, 2차 p=0.99). 임령 자체가 다양성 수준을 정하지 않음 — DHI 생산성 혹형과 대조되는 정직한 null.

**(ii) 임령이 RS–다양성 관계를 조절 — 11/19 유의, 전부 음(−):**
| 관계 | 상호작용 | p |
|---|---|---|
| Hill q1 ~ LAI × age | −0.214 | 0.002 |
| Hill q1 ~ VCI × age | −0.168 | 0.003 |
| Hill q1 ~ SAVI × age | −0.179 | 0.005 |
| Hill q2 ~ EVI × age | −0.178 | 0.007 |
| Turnover ~ EVI × age | −0.146 | 0.012 |
- **핵심**: RS(분광 SAVI/EVI/NDVI + 구조 LAI/VCI/FHD)–다양성 커플링은 **젊은 임분에서 강하고 성숙할수록 약해짐**. 성숙 폐쇄 캐노피에서 분광·구조 신호가 포화→종다양성과 디커플링. **님 생장 가설의 핵심 귀결.**

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
  - Canopy_Ht_trend는 bivariate 음→통제 후 양 **부호 반전**(Simpson's paradox: 사이트 내에선 노령일수록 수고 trend 양).

### 30m 재실행 (robustness, 검열 감안)
30m(≤40 검열)으로 (B) 재실행 시 조절효과 패턴이 GAMI와 갈림(일부 반대 부호). **모순 아님** — 30m은 76% 포화라 사실상 "최근교란된 젊은 25%" 대비이므로 진짜 임령이 아닌 recency 축을 측정. 따라서 **(B)의 유효 결론은 GAMI(진짜 임령) 기준**이며, 30m은 젊은층 검증·recency proxy로만 사용. (Obj2 severity/recency 결과와 정합)

### 종합 (한 문장)
> 임령은 다양성 *수준*을 직접 정하진 않지만, **RS–종다양성 관계의 강도(젊을수록 강)와 구조 복잡화의 양상(젊음=축적, 노령=갭동태)을 좌우**한다. 즉 원격탐사로 다양성을 추론할 때 **천이단계가 조절자**다.

그림 `N01_stand_age.png`. 결과 `stand_age_{models,moderation,trendage}_gami.csv`.

## 재현
```
python scripts/obj23/130_trend_increment.py          # (A)
python scripts/obj23/141_extract_age_30m.py          # 30m 타일 추출 (다운로드 후)
python scripts/obj23/140_stand_age.py                # (B) 혼합모델, 존재하는 age 소스 자동 병합
python scripts/obj23/142_stand_age_figures.py        # N02 교차검증 + N03 임령~구조변화 산점도
```
