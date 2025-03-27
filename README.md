# Description
- **Title**: Developed simulation and R-based integrated analysis tools to detect patterns of variability differences in data   
**`데이터 내 변동성 차이 패턴 검출을 위한 시뮬레이션 및 R 기반 통합 분석 도구 개발`**
  
- **Tools**:&nbsp;&nbsp;
![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=R&logoColor=white)
![Illustrator](https://img.shields.io/badge/Illustrator-FF9A00?style=flat-square&logo=adobeillustrator&logoColor=white)
- **Skills**: Simulation study, Statistical modeling, Method modification, Parallel processing, Data preprocessing, Visualization

### 📄 Publication Info
- **Master’s Thesis**: Software development and application to identify differentially dispersed genes based on transcriptome data (2025)
- **Institution**: Korea University, Graduate School of Computer and Information Science
<br><br>
# Summary
### 🤔 Problem

암과 같이 예측이 어려운 질병은 유전자 발현이 불규칙하게 나타나는 반면, 감기처럼 일시적인 질환에서는 발현이 비교적 일관됨. 이처럼 질병 상태에 따른 유전자 발현의 변동성 차이는 생물학적 해석의 중요한 단서가 될 수 있지만, 기존 연구는 특정 분석 방법에 의존할 뿐 체계적인 연구는 여전히 부족함

### 🔍 Approach & Workflow

이 문제를 해결하기 위해, 아래와 같은 분석 프로세스를 설계하였음

**A. 실제 생물학적 데이터가 가질 수 있는 다양한 상황을 반영하는 가상 시뮬레이션 데이터 생성**

- WEHI 연구기관의 통계 그룹에서 good-turing 알고리즘을 기반으로 추정한 가상 데이터 생성 방식을 활용함
 
- 본 연구 목적에 맞게 데이터 내 **9가지 변인을 조절**하여 (ex. 샘플 수, 그룹 간 샘플 수 불균형, 정규화 방식 등) 10,000개의 차원을 가지는 **가상 데이터를 약 62,000개 생성**함 


**B. 패턴 추정 가능한 통계 방법 조사 (49가지) ➔ 적용 가능성 검사 및 필터링 (23+1가지) ➔ 일부 통계 방법은 분석 가능한 형태로 수정**



**C. A에서 생성한 모든 가상 데이터에 24가지 분석 방법 모두 적용 후 벤치마크 비교/분석**

데이터 특성에 따라 최적의 방법 도출 


**D. 실제 생물학적 데이터에 적용 및 실증**


### 📈 Conclusion
위 과정을 토대로 실제 데이터 input 시 데이터를 분석에 용이하게 전처리, 패턴 검출 분석 방법 적용, 결과를 제공해주는 통합 분석 도구를 개발함 
<br><br>
# File Structure

✅ Data .

### 📁 Data 
Supplementary Table들 넣기 = 시뮬레이션 결과


### 📑 Scripts
R 디렉토리 = 패키지 개발내 들어가는 모든 것을 모아둔 R 전부 업로드




















