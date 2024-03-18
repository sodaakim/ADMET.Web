# ADMET 분석 서비스용 웹서버

![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/89abc3cf-dcbe-4eb6-a134-cc8942525c02)

ADMET 모델은 분자의 물리화학적 성질, 독성 등을 분석할 수 있어 신약 개발 과정에서 중요한 도구로 활용되며, 편리한 연구 증진을 위해 Flask와 MySQL로 ADMET 모델을 제공하는 웹서비스를 개발하고 있습니다.

SMILES로 표기된 화합물을 입력하면 Property 86개에 대한 분석 결과를 제공합니다.

연구자들의 경험을 반영한 사용자 친화적 웹 애플리케이션이며, 본 프로젝트는 확장성과 유지보수를 고려한 컴포넌트 기반 설계와 서버의 처리 한계를 고려한 대용량 데이터 관리 기법을 특징으로 합니다.

This tool supports researchers in identifying and optimizing potential risk factors by predicting the absorption, distribution, metabolism, excretion, and toxicity (ADMET) profiles of drug candidate substances. It offers a user-friendly interface, enabling scientists to easily predict the bioavailability, safety, and efficacy of drugs.

## 기능 사항
- **SMILES 입력**

단일 SMILES 입력(Analyze)과 다중 SMILES 입력(Screening)이 가능합니다.

![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/01116ebc-f084-49fa-825d-159db74e155d)

- **Property 분석**

86개의 ADMET Property에 대해 분석합니다.

ADMET 분석을 위한 모델은, 본 연구소에서 개발한 모델 뿐만 아니라 FDADMET와 OPERA의 분석 모델이 포함되어있습니다.

![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/be3f2475-f177-4d1e-97ac-e903a9141293)

- **결과 다운로드**

PDF와 CSV로 결과를 다운로드 할 수 있습니다.

다운로드 옵션에서 Property를 선택해 원하는 사항만 다운로드 할 수 있습니다.

![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/2e489392-ebe5-4808-ad27-dc38284d942a)

![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/dce705eb-962b-4fcc-a1b9-5c8619f639c8)
![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/0eca17a0-3609-4c85-840d-6ac2537842d1)


## 개발 과정 하이라이트
- **사용자 친화적 디자인** : 연구자의 경험을 반영하여 최적화된 UI/UX.
- **모듈화** : Flask Blueprint를 활용한 클린 아키텍처 구현.
- **컴포넌트 기반 설계** : 재사용 가능한 UI 컴포넌트를 통한 높은 유지보수성.
- **대용량 데이터 처리** : 서버의 한계를 고려한 효율적인 데이터 관리와 분석.

## 공개 사항
개발중인 웹 서비스의 보안을 위해 일부 기능만을 공개하고 있습니다.



## 시작하기
이 섹션에서는 프로젝트를 로컬 환경에서 실행하기 위한 단계를 설명합니다.

### 전제 조건
- Python 3.8+
- requirements.txt 에 포함된 실행 요구조건을 확인해주세요.

### 설치
```bash
git clone https://github.com/sodaakim/ADMET.Web.git
cd ADMET.Web
pip install -r requirements.txt
```

### 사용 방법
```bash
flask run
```

## 연구소 정보
중앙대학교 융합공학부 바이오메디컬공학 생명정보학 연구실

Department of Biomedical Engineering

School of integrative engineering Chung-Ang University

84 Heukseok-ro, Dongjak-gu, Seoul, Republic of Korea

### Contact Us

- **Homepage** : http://ssbio.cau.ac.kr


### 개발자 정보
- **이름** : 김소희

- **이메일** : sohuikim2020@naver.com
