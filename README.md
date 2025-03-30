# Description
- **Title**: Develop simulation and R-based integrated analysis tool to detect pattern of variability differences in data   
**`데이터 내 변동성 차이 패턴 검출을 위한 시뮬레이션 및 R 기반 통합 분석 도구 개발`**
  
- **Tools**:&nbsp;&nbsp;
![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=R&logoColor=white)
![Illustrator](https://img.shields.io/badge/Illustrator-FF9A00?style=flat-square&logo=adobeillustrator&logoColor=white)
- **Skills**: Statistical simulation, Statistical modeling, Parallel processing, Data preprocessing, Visualization

### 📄 Publication Info
- **Master’s Thesis**: Software development and application to identify differentially dispersed genes based on transcriptome data (2025)
- **Institution**: Korea University, Graduate School of Computer and Information Science
<br><br>
# Summary
### 🤔 Problem

건강한 사람들의 유전자 발현은 일반적으로 일정한 범위 내에서 유지되지만, 암과 같이 예측이 어려운 질병에서는 환자마다 발현이 크게 달라지는 경우가 많음. 이처럼 질병 상태에 따라 발생하는 **발현의 불안정성 = 변동성(분산)의 차이**는 생물학적 해석에 중요한 단서가 될 수 있지만, 기존 연구는 특정 분석 방법에 의존할 뿐 체계적인 연구가 부족함 

<img src="figures/Target pattern.png" width="900"/>

### 🔍 Approach & Workflow

이 문제를 해결하기 위해, 아래와 같은 분석 프로세스를 설계하였음

<br>

1️⃣ **실제 생물학적 데이터가 가질 수 있는 다양한 상황을 반영하는 가상 시뮬레이션 데이터 생성**

- WEHI 연구기관의 통계 그룹에서 good-turing 빈도 추정 방식을 기반으로 구축한 가상 생물학 데이터 생성 방식을 활용함
 
- 본 연구 목적에 맞게 데이터 내 **9가지 조건을 조절**하여 (ex. 샘플 수, 그룹 간 샘플 수 불균형, 정규화 방식 등) 10,000개의 차원을 가지는 **가상 데이터를 약 62,000개 생성**함 

<br>

2️⃣ **변동성 차이 추정이 가능한 통계 방법 정립**

- 전통적 통계 방법, 선형 회귀 기반 방법, 특정 분야에서 개발된 방법 등 **변동성 차이를 추정할 가능성 있는 다양한 방법 조사 (49가지)**
  
- 소규모 테스트 데이터를 통해 적용 가능성 검사 및 성능이 낮거나 비효율적인 방법 **사전 필터링 (23+1가지)**
  
- 일부 통계 방법에 대해 **분석 가능한 형태로 방법론 코드 수정**
<br>
<img src="figures/Method used.png" width="800"/>

<br>

3️⃣ **모든 가상 데이터에 24가지 분석 방법 적용 후 성능 비교/분석**

그 결과,

- 방법론은 크게 두 그룹으로 분류됨 (통계적 검정력과 제 1종 오류가 모두 높은 과대추정 성향 또는 모두 상대적으로 낮은 보수적 성향)

- Zero-inflated 문제를 분포 수준에서 모델링하여 보정하는 MDSeq은, 일정 샘플 수 이하 조건에서 가장 우수한 성능을 보임

- 선형 회귀 모델의 오차항에 보조 회귀를 적용하는 White’s 검정은 어떤 조건에서도 제1종 오류를 가장 안정적으로 억제하며, 샘플 수가 증가할수록 다른 방법에 준하는 검정력을 보임

- 이처럼 **각 9가지 서로 다른 조건에서 가장 성능이 안정적인 최적의 방법을 도출**함
<br>
<img src="figures/Simulation result example 1.jpg" width="800"/>


### 📈 Conclusion

✅ **R 기반 통합 분석 도구 개발**

아래와 같은 기능을 제공하는 도구를 개발함
- 실제 데이터의 특성을 파악하여 분석에 용이하게 자동 전처리
- **병렬 처리 기반 최적의 변동성 차이 분석** 수행
- 분석 결과 통합 정리 및 시각화 보조를 구현 결과 구조화
<br>

✅ **실제 생물학적 데이터에서 도구의 필요성 실증**
- 미국 COPD 환자 코호트의 생물학 데이터를 바탕으로 흡연자와 금연자 간 최적의 변동성 차이 분석 방법을 적용하여 본 연구의 필요성을 실증함
     
- 추가적으로 변동성 차이 패턴에 대한 **생물학적 해석을 위한 알고리즘을 제안 및 활용**하여 새로운 인사이트를 도출할 수 있음을 증명함

<img src="figures/Substantiation.jpg" width="600"/>



<br><br>
# File Structure

모든 시뮬레이션 벤치마크 결과 및 툴 사용을 위한 R 스크립트가 업로드되었으며, 시각화용 스크립트는 업로드 하지 않음

### 📁 Result 
**Supplementary Table S1.xlsx**: &nbsp;모든 시뮬레이션 디자인 별 조절한 파라미터 정보

**Supplementary Table S2.xlsx**: &nbsp;**`기본 파라미터 2개를 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S3.xlsx**: &nbsp;**`기술적 변동계수를 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S4.xlsx**: &nbsp;**`임의의 노이즈를 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S5.xlsx**: &nbsp;**`라이브러리 크기의 불균형을 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S6.xlsx**: &nbsp;**`다양한 정규화 방식을 적용한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S7.xlsx**: &nbsp;**`그룹 간 샘플 크기의 불균형을 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S8.xlsx**: &nbsp;**`상대적 발현 수준을 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S9.xlsx**: &nbsp;**`0의 비율을 조절한 데이터`** 에 대한 24가지 통계 검정의 성능 평가 결과

**Supplementary Table S10.xlsx**: &nbsp;**실제 미국 COPD 환자 코호트**에서 흡연자와 금연자 간 유전자 발현 **변동성 차이 분석** 결과 

**Supplementary Table S11.xlsx**: &nbsp;흡연자와 금연자 간 발현 변동성 차이가 나는 유전자들의 생물학적 기능에 대한 **Fisher's exact test** 결과


### 📑 R scripts

**DvarSeq.R**: &nbsp; 데이터 내 변동성 차이 패턴을 검출하는 통합 분석 도구를 구현한 스크립트

**fileLoad.R**: &nbsp; 분석 도구 실행에 필요한 패키지 및 함수들을 로드하는 스크립트

**simGen.R**: &nbsp; 다양한 조건을 반영하여 약 62,000개의 가상 데이터를 생성하는 스크립트

**simDvarSeq.R**: &nbsp; 생성된 가상 데이터 전체에 DvarSeq 함수를 병렬로 적용하는 스크립트

**simControl.R**: &nbsp; 9가지 조건 설정부터 데이터 생성, 변동성 분석까지 전체 과정을 제어하는 메인 스크립트

**Other R scripts**: &nbsp; 그 외 모든 파일은 분석 도구에 포함된 다양한 분석 방법 및 관련 보조 함수들을 구현한 스크립트로, 외부 소스에서 복사하거나 일부 수정하여 사용하였으며, 패키지 의존성을 줄이기 위해 직접 코드를 포함함


### 🔓 License

This repository includes original code and multiple functions adapted from various sources under GPL-3 and compatible licenses (e.g., skedastic). See individual R files for function-level attribution.
















