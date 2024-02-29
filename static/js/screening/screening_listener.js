document.addEventListener('DOMContentLoaded', function() {
    let smilesArray;

    // "Send" 버튼 이벤트 리스너
    document.querySelector('.button.send').addEventListener('click', function() {
        let smilesInput = document.getElementById('smilesInput').value;
        if (!smilesInput.trim()) {
            alert("Please input SMILES strings.");
            return;
        }

        // 줄바꿈으로 구분하여 SMILES 배열 생성
        smilesArray = smilesInput.split('\n').filter(smile => smile.trim() !== '');

        // 유효한 SMILES 문자열의 개수를 출력
        document.getElementById('Valid_molecules').textContent = smilesArray.length;
    });

    // "Result" 버튼 이벤트 리스너
    document.querySelector('.button.result').addEventListener('click', function() {
        if (!smilesArray || smilesArray.length === 0) {
            alert("Please enter SMILES strings and click Send first.");
            return;
        }

        // 분석 시작 시 로딩 오버레이 표시
        document.getElementById('loadingOverlay').style.display = 'flex';

        // 서버에 분석을 요청
        fetch('/screening', { // '/screening' 엔드포인트로 POST 요청 보냄
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ smiles_list: smilesArray })
        })
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            // 응답이 성공적이면, 사용자를 결과 페이지로 리디렉션
            window.location.href = '/result'; // 쿼리 파라미터 '?page=1'은 필요에 따라 추가
        })
        .catch(error => {
            console.error('Error:', error);
            alert('There was an error processing your request.');
        })
        .finally(() => {
            // 로딩 오버레이 숨김
            document.getElementById('loadingOverlay').style.display = 'none';
        });
    });

    /*
    document.querySelector('.button.result').addEventListener('click', function() {
        if (!smilesArray || smilesArray.length === 0) {
            alert("Please enter SMILES strings and click Send first.");
            return;
        }

        // 분석 시작 시 로딩 오버레이 표시
        document.getElementById('loadingOverlay').style.display = 'flex';

        // 서버에 분석을 요청
        let fetchPromise = fetch('/screening', { // '/screening' 엔드포인트로 POST 요청 보냄
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ smiles_list: smilesArray })
        });
        let timeoutPromise = new Promise(resolve => setTimeout(resolve, 3000));

        Promise.all([fetchPromise, timeoutPromise]).then(async ([response]) => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            // 응답이 성공적이면, 사용자를 결과 페이지로 리디렉션
            window.location.href = '/result'; // 쿼리 파라미터 '?page=1'은 필요에 따라 추가
        })
        .catch(error => {
            console.error('Error:', error);
            alert('There was an error processing your request.');
        })
        .finally(() => {
            // 로딩 오버레이 숨김
            document.getElementById('loadingOverlay').style.display = 'none';
        });
    });*/

    document.addEventListener('click', function(event) {
        if (event.target.classList.contains('detail')) {
            // Detail 버튼이 클릭된 경우
            const cardElement = event.target.closest('.card');
            if (cardElement) {
                const smiles = cardElement.querySelector('.card-contents h4:nth-of-type(2)').textContent;
                alert(`Detail for SMILES: ${smiles}`);
                // 여기서 상세 정보를 표시하는 로직을 구현할 수 있음
                // 예: 상세 페이지로 이동하거나, 모달 창을 표시하는 등
            }
        }
    });
});
