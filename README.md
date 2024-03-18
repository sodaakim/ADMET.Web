# ADMET 분석 서비스용 웹서버
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 236.88 115.44">
  <defs>
    <style>
    </style>
    <filter id="drop-shadow-1" filterUnits="userSpaceOnUse">
      <feOffset dx="0" dy="0"/>
      <feGaussianBlur result="blur" stdDeviation="1.42"/>
      <feFlood flood-color="#231815" flood-opacity=".4"/>
      <feComposite in2="blur" operator="in"/>
      <feComposite in="SourceGraphic"/>
    </filter>
  </defs>
</svg>![image](https://github.com/sodaakim/ADMET.Web/assets/83997634/2164aab4-cdfa-4edd-baa0-cc64d69e3cef)


ADMET 모델은 분자의 물리화학적 성질, 독성 등을 분석할 수 있어 신약 개발 과정에서 중요한 도구로 활용되며, 편리한 연구를 위해 Flask와 MySQL로 웹서비스를 개발하고 있습니다.

SMILES로 표기된 화합물을 입력하면 Property 86개에 대한 분석 결과를 제공합니다.

연구자들의 경험을 반영한 사용자 친화적 웹 애플리케이션이며, 본 프로젝트는 확장성과 유지보수를 고려한 컴포넌트 기반 설계와 서버의 처리 한계를 고려한 대용량 데이터 관리 기법을 특징으로 합니다.

This tool supports researchers in identifying and optimizing potential risk factors by predicting the absorption, distribution, metabolism, excretion, and toxicity (ADMET) profiles of drug candidate substances. It offers a user-friendly interface, enabling scientists to easily predict the bioavailability, safety, and efficacy of drugs.

## 기능 사항
- **SMILES 입력**

- 


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

- **Email** : blisszen@cau.ac.kr

### 개발자 정보
- **이름** : 김소희

- **이메일** : sohuikim2020@naver.com
