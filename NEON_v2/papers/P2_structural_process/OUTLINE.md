# P2 아웃라인 — 프로세스: 다중시기 LiDAR로 본 산림 구조 발달

## 제목(안)
- (주) "Multitemporal LiDAR reveals divergent successional and disturbance signatures of forest structural change"
- (부/국문) 다중시기 LiDAR로 포착한 산림 구조 발달 — 천이와 교란의 상반된 서명

## 타깃 저널
Journal of Ecology / Forest Ecology and Management / Remote Sensing of Environment

## 초록 (1문단 초안)
다중시기 LiDAR로 캐노피 구조 변화율(ΔLiDAR/연간 trend)을 측정해 두 동인 — **천이(임분연령)와 교란(press/pulse)** — 을 규명했다. 천이 축에서 젊은 임분은 수직 복잡도·엽면적을 빠르게 축적하고(VCI/LAI trend −0.25/−0.22 vs 임령), 노령림은 갭이 열리는 gap-phase 동태로 전환되었다(Deep_Gap trend +0.28). 교란 축에서 press(곤충 고사, ΔLAI −0.18)와 pulse(산불·풍해 재생, ΔLAI +0.22)는 상반된 구조 서명을 보였다. 평균회귀를 엄밀 통제하자 초기 복잡도의 완충은 **Deep_Gap·Gini에서만** 진짜였다. 원격탐사가 산림 구조 발달의 천이·교란 동태를 분해함을 보인다.

## 핵심 주장(기여)
1. **천이 구조 궤적**: aggradation(복잡화 축적) → old-growth(갭-단계). 임령으로 예측(부분회귀 생존 4지표).
2. **press vs pulse**: 다중시기 LiDAR가 교란 체제 구분(Choi 2023 지지).
3. **복잡도 완충 = RTM 통제 후 Deep_Gap·Gini만**: naive 완충의 대부분은 평균회귀 아티팩트.

## Figure 순서 (main)
1. **N04** — 순수 임령효과(부분회귀): VCI/LAI −, Deep_Gap/Canopy_Ht + ★
2. **N03** — 임령 vs 구조변화 산점도(site색, plateau)
3. **N06** — Canopy_Ht Simpson 분해(between vs within)
4. **J01** — press vs pulse 구조 서명
5. **J05** — RTM 통제 후 생존 지표(Deep_Gap·Gini)
- Supp: N01·N02(임령 교차검증 30m 포화)·N05·N07·N08, J02~J04·J06

## 구성
- Intro: RS로 산림 구조 발달 읽기 / 천이 + 교란 두 동인 / stand development 이론
- Methods: 다중시기 LiDAR 지표, GAMI 100m 임령(+30m 교차검증), ΔLiDAR/trend, RTM 통제(impact×complexity), Choi2023 press/pulse
- Results: (1)천이 구조변화 (2)press/pulse (3)복잡도 완충 RTM
- Discussion: aggradation→gap-phase, RTM 방법론 교훈, RS가 구조 동태 포착

## 한계 (정직)
- **임령 사이트단위(ICC=0.75)**, GAMI 100m 아티팩트(ABBY 내부SD 70년) → 유효 N≈19, 사이트 수준 해석
- FHD_trend 등은 통제 시 소멸(생물지리) → 순수 임령효과는 VCI/LAI/Deep_Gap/Canopy_Ht 4개
- LiDAR 관측창 ~8–10년(단기), 30m 임령은 recency proxy(≤40 포화)
