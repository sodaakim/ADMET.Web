document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('exampleButton').addEventListener('click', function() {
        document.getElementById('smilesInput').value = 'C1=CC=CC=C1O'; // Example SMILES
    });
});


document.addEventListener('DOMContentLoaded', function() {
    /*
    function hideLoadingOverlay() {
        document.getElementById('loadingOverlay').style.display = 'none';
    }
    // 페이지 로드와 포커스를 얻을 때 로딩 오버레이 숨기기
    hideLoadingOverlay();
    window.addEventListener('focus', hideLoadingOverlay);*/

    document.querySelector('.button.analyze').addEventListener('click', function() {
        let smiles = document.getElementById('smilesInput').value;

        if (!smiles) {
            alert("Please enter a SMILES string.");
            return;
        }

        // 분석 시작 시 로딩 오버레이 표시
        document.getElementById('loadingOverlay').style.display = 'flex';
        // 서버 요청을 시작하는 Promise
        let fetchPromise = fetch('/analyze', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ smiles: smiles })
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
            if (data.image_url) {
                // 이미지 URL과 SMILES 문자열을 /evaluation 페이지로 리디렉션
                window.location.href = `/evaluation?image_url=${encodeURIComponent(data.image_url)}&smiles=${encodeURIComponent(smiles)}`;
                //window.location.href = `/evaluation?smiles=${encodeURIComponent(smiles)}`;
            } else {
                alert('There was an error processing the SMILES string.');
            }
        }).catch(error => {
            console.error('Error:', error);
            alert('There was an error processing your request.');
        }).finally(() => {
            // 모든 작업이 끝난 후 로딩 오버레이 숨김
            document.getElementById('loadingOverlay').style.display = 'none';
        });
    });
});


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
});