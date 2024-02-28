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
            alert("Please enter a SMILES string.");
            return;
        }

        // 분석 시작 시 로딩 오버레이 표시
        document.getElementById('loadingOverlay').style.display = 'flex';

        let fetchPromise = fetch('/result', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ smiles_list: smilesArray })
        });

        // 최소 3초 대기를 보장하는 Promise
        let timeoutPromise = new Promise(resolve => setTimeout(resolve, 3000));

        // 서버 요청과 3초 대기 중 늦게 완료되는 것을 기다림
        Promise.all([fetchPromise, timeoutPromise]).then(async ([response]) => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.json();
        }).then(data => {
            // 서버로부터의 응답에 따라 적절한 리디렉션 수행
            if (data && data.length > 0 && data[0].image_url) {
                // 예시: 첫 번째 SMILES 문자열에 대한 결과 페이지로 리디렉션
                window.location.href = `/evaluation?image_url=${encodeURIComponent(data[0].image_url)}&smiles=${encodeURIComponent(smilesArray[0])}`;
            } else {
                alert('There was an error processing the SMILES string.');
            }
        }).catch(error => {
            console.error('Error:', error);
            alert('There was an error processing your request.');
        }).finally(() => {
            // 로딩 오버레이 숨김
            document.getElementById('loadingOverlay').style.display = 'none';
        });
    });
});


/*
document.addEventListener('DOMContentLoaded', function() {
    // URL에서 쿼리 매개변수를 파싱하는 함수
    function getQueryVariable(variable) {
        var query = window.location.search.substring(1);
        var vars = query.split('&');
        for (var i = 0; i < vars.length; i++) {
            var pair = vars[i].split('=');
            if (decodeURIComponent(pair[0]) == variable) {
                return decodeURIComponent(pair[1]);
            }
        }
        return false;
    }

    // SMILES 문자열을 가져와서 페이지에 설정
    var smiles = getQueryVariable('smiles');
    if (smiles) {
        document.querySelector('.molecule-smiles .smiles').textContent = smiles;
    } else {
        console.log("SMILES string not found in URL");
    }

    // 이미지 URL을 가져와서 페이지에 설정
    var image_url = getQueryVariable('image_url');
    if (image_url) {
        document.querySelector('.molecule-image .molecule').src = image_url;
    } else {
        console.log("Image URL not found in URL");
    }
});*/