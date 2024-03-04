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
            window.location.href = `/evaluation?smiles=${encodeURIComponent(smiles)}`;
        }).catch(error => {
            console.error('Error:', error);
            alert('There was an error processing your request.');
        }).finally(() => {
            // 모든 작업이 끝난 후 로딩 오버레이 숨김
            document.getElementById('loadingOverlay').style.display = 'none';
        });
    });
});
