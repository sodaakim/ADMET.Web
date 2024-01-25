body {
    font-family: 'Arial', sans-serif;
    margin: 0;
    padding: 0;
    box-sizing: border-box;
    color: #333;
    background: #f4f4f4; /* 페이지 배경색 변경 */
}

header {
    background-color: #005f73; /* 진한 녹색으로 변경 */
    color: white;
    padding: 15px 0;
    box-shadow: 0 2px 5px rgba(0, 0, 0, 0.2); /* 상단에 그림자 효과 추가 */
    position: sticky; /* 헤더를 상단에 고정 */
    top: 0;
    z-index: 1000;
}

header .logo a {
    font-size: 26px; /* 로고 크기 증가 */
    text-decoration: none;
    font-weight: bold;
    color: white;
    transition: color 0.3s ease; /* 색상 변화 애니메이션 */
    padding-left: 25px;
}

header .logo a:hover {
    color: #ade8f4; /* 마우스 오버 시 연한 녹색으로 변경 */
}

nav ul {
    list-style-type: none;
    margin: 0;
    padding: 0;
    display: flex;
    justify-content: center;
    align-items: center;
}

nav ul li {
    margin: 0 15px; /* 메뉴 항목 간격 조정 */
    position: relative;
}

nav ul li a {
    text-decoration: none;
    color: white;
    padding: 10px 15px;
    display: inline-block;
    transition: all 0.3s ease; /* 전체 애니메이션 효과 적용 */
}

nav ul li a:hover,
nav ul li a:focus { /* 포커스 추가로 접근성 향상 */
    background-color: #0077b6; /* 호버 시 배경색 변경 */
    border-radius: 4px; /* 둥근 모서리 추가 */
    transform: scale(1.05); /* 약간 확대 */
}

/* 드롭다운 메뉴 스타일 */
.dropdown-content {
    display: none;
    position: absolute;
    left: 50%;
    transform: translateX(-50%); /* 센터 정렬 */
    background-color: white;
    min-width: 200px; /* 드롭다운 너비 조정 */
    box-shadow: 0 8px 16px rgba(0,0,0,0.1); /* 드롭다운 그림자 추가 */
    border-radius: 4px; /* 드롭다운 둥근 모서리 추가 */
    z-index: 1;
}

.dropdown-content a {
    color: #005f73; /* 드롭다운 메뉴 아이템 색상 변경 */
    padding: 12px 16px;
    text-decoration: none;
    display: block;
}

.dropdown-content a:hover {
    background-color: #0a9396; /* 드롭다운 메뉴 아이템 호버색 변경 */
    color: white;
}

.dropdown:hover .dropdown-content {
    display: block;
}

/* 반응형 메뉴 디자인 */
@media (max-width: 768px) {
    header {
        padding: 10px 0;
    }

    header .logo a {
        font-size: 24px;
    }

    nav ul {
        flex-direction: column;
    }

    nav ul li {
        margin: 10px 0;
    }

    .dropdown-content {
        left: 0;
        transform: none;
    }
}