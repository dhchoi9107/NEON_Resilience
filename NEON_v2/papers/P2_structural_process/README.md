# P2 — 프로세스: 다중시기 LiDAR로 본 산림 구조 발달 (천이 + 교란)

## 주제
**구조 다양성의 변화(ΔLiDAR / 연간 trend)**를 두 동인으로 설명: **천이(임분연령)** 와 **교란(press/pulse)**. 원격탐사(다중시기 LiDAR)가 산림 구조 발달 동태를 포착함을 규명.

## 핵심 결과
1. **천이 구조변화** (`N03·N04`, GAMI 임령): 젊을수록 **VCI·LAI trend 큼**(수직복잡도·엽면적 빠른 축적, β=−0.25/−0.22), 노령림 **Deep_Gap trend 양**(갭-단계 동태, +0.28), Canopy_Ht 부호반전(Simpson, `N06`). = 교과서적 stand development(aggradation→gap-phase).
   - **정직한 한계**: 임령 ICC=0.75(사이트 단위), FHD_trend는 통제 시 소멸(교란). 순수 임령효과는 부분회귀(N04)로 VCI/LAI/Deep_Gap/Canopy_Ht 4개만 생존. 유효 N≈19.
2. **교란 구조변화** (`J01`, Choi 2023): press(곤충)=ΔLAI −0.18(고사) vs pulse(산불/풍해)=ΔLAI +0.22(재생). 다중시기 LiDAR가 교란 체제 구분.
3. **복잡도 완충 (RTM 통제)** (`J04·J05`): naive 완충은 대부분 평균회귀. 엄밀 통제 후 **Deep_Gap(완충)·Gini(증폭)만 진짜**, 수고·LAI·FHD는 RTM.

## 그림 (figures/)
- **N01** 임령~복잡화 / **N02** GAMI vs 30m 교차검증(30m 40년 포화)
- **N03** 임령~구조변화 산점도(site색) / **N04** 부분회귀(순수 임령효과) ★
- **N05** bivar vs 부분회귀(도메인 묶음) / **N06** Canopy_Ht Simpson 분해
- **N07·N08** 임령×구조변화 점크기=다양성(raw/partial)
- **J01** press/pulse / **J02·J03** 복잡도 완충·클래스 / **J04·J05·J06** RTM 통제·생존·계산

## 결과 (results/)
stand_age_trendage_{gami,30m}, choi2023, complexity_rtm, stand_age_crossval

## 데이터
`../../data/plot_stand_age_{gami,30m}.csv` (GAMI 100m = 유효 임령; 30m = recency proxy ≤40 포화)

## 관련 스크립트 (../../scripts/obj23/)
`140·141`(임령 모델·30m 추출)·`142`(N02~N06)·`143·144`(N07·N08)·`90·91·92`(Choi2023·RTM)
