function updatePagination() {
    const paginationElement = document.getElementById('pagination');
    paginationElement.innerHTML = ''; // 기존 페이지네이션 버튼을 초기화

    // 페이지네이션 버튼 생성
    for (let i = 1; i <= totalPages; i++) {
        // 현재 페이지 주변의 페이지 번호만 보여주기 (예: 현재 페이지 +-2)
        if (i === currentPage || i === currentPage - 1 || i === currentPage + 1 || i === 1 || i === totalPages) {
            const pageButton = document.createElement('button');
            pageButton.textContent = i;
            pageButton.className = 'button page-button';
            pageButton.disabled = i === currentPage; // 현재 페이지 버튼 비활성화
            pageButton.addEventListener('click', function() {
                currentPage = i;
                showPage(currentPage);
            });
            paginationElement.appendChild(pageButton);

            // ...와 같은 구분자 추가
            if (i === currentPage - 2 && i > 2) {
                const dots = document.createElement('span');
                dots.textContent = '...';
                paginationElement.insertBefore(dots, pageButton);
            } else if (i === currentPage + 2 && i < totalPages - 1) {
                const dots = document.createElement('span');
                dots.textContent = '...';
                paginationElement.appendChild(dots);
            }
        }
    }

    // 'Previous'와 'Next' 버튼의 상태 업데이트
    updatePaginationButtons();
}

// 'showPage' 함수 내에서 'updatePagination' 함수 호출
function showPage(page) {
    // (카드 표시 로직...)

    updatePagination(); // 페이지네이션 업데이트
}

// 페이지 로드 시 초기 페이지 설정 및 페이지네이션 구성
document.addEventListener('DOMContentLoaded', function() {
    showPage(currentPage);
});
