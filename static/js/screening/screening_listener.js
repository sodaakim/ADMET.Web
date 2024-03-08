document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('exampleButton').addEventListener('click', function() {
        document.getElementById('smilesInput').value =
                                                                'OOS(=O)(=O)O'
                                                                +'\n'+'NOS(=O)(=O)O'
                                                                +'\n'+'COCC(CO)O'
                                                                +'\n'+'COC(=O)C=C.C=C'
                                                                +'\n'+'CN=C1CN(C(=C2C=C(C=CC2=N1)Cl)C3=CC=CC=C3)O'
                                                                +'\n'+'CN(C)CCN(C)C'
                                                                +'\n'+'CN(C)C1=CC=CC(=C1)C2=CC(=O)C3=C(N2)C=CC(=C3)OC'
                                                                +'\n'+'CCSCCO'
                                                                +'\n'+'CCOC(=O)Cl'
                                                                +'\n'+'CCCN1C(=O)C2=CC=CC=C2N=C1SCC3=NC4=C(C(=C(S4)C)C)C(=O)N3'
                                                                +'\n'+'CCCCNC1=CC=C(C=C1)C(=O)O'
                                                                +'\n'+'CCCC[C@@H](C(C)C)OCC'
                                                                +'\n'+'CCC(CC)C(=O)OC(C)C'
                                                                +'\n'+'CCC(C)CCC(=O)OCC'
                                                                +'\n'+'CC(F)(Cl)Cl'
                                                                +'\n'+'CC(C)CC(=O)OCC(C)C'
                                                                +'\n'+'CC(C)(C)OCC(C)(C=C)O'
                                                                +'\n'+'CC(=O)OCCOC'
                                                                +'\n'+'CC(=O)OC=C.C=C'
                                                                +'\n'+'CC(=CCCC(OC)OC)C'
                                                                +'\n'+'C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2'
                                                                +'\n'+'C1=CC2=C3C(=C1)C=C4C=CC5=C4C3=C(C=C2)C=C5'
                                                                +'\n'+'C1=CC=C2C(=C1)N=C(O2)SCCS(=O)(=O)C3=NC4=CC=CC=C4S3'
                                                                +'\n'+'C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C(=C(C=C3)O)O'
                                                                +'\n'+'C1=CC=C(C=C1)COC2=CC3=CC=CC=C3C=C2'
                                                                +'\n'+'C1=CC=C(C=C1)CN2C=C(C3=CC=CC=C32)CO'
                                                                +'\n'+'C[N+](C)(C)CCO'
                                                                +'\n'+'C(CSCC(CS)SCCS)S'
                                                                +'\n'+'C(C(Cl)(Cl)Cl)(Cl)(Cl)Cl'
                                                                +'\n'+'C(C(CBr)Br)Cl'+'\n'; // Example SMILES
    });
});


let smilesArray;

document.addEventListener('DOMContentLoaded', function() {

    /*
    // "Send" 버튼 이벤트 리스너
    document.querySelector('.button.send').addEventListener('click', function() {
        let smilesInput = document.getElementById('smilesInput').value;
        if (!smilesInput.trim()) {
            alert("Please input SMILES strings.");
            return;
        }

        smilesArray = smilesInput.split('\n').filter(smile => smile.trim() !== '');
        localStorage.setItem('smilesArray', JSON.stringify(smilesArray));

        // 유효한 SMILES 문자열의 개수를 출력
        document.getElementById('Valid_molecules').textContent = smilesArray.length;
    }); */

    // "Result" 버튼 이벤트 리스너
    document.querySelector('.button.result').addEventListener('click', function() {

        let smilesInput = document.getElementById('smilesInput').value;
        if (!smilesInput.trim()) {
            alert("Please input SMILES strings.");
            return;
        }

        smilesArray = smilesInput.split('\n').filter(smile => smile.trim() !== '');
        localStorage.setItem('smilesArray', JSON.stringify(smilesArray));

        // 유효한 SMILES 문자열의 개수를 출력
        // document.getElementById('Valid_molecules').textContent = smilesArray.length;


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
            window.location.href = '/result';
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
});


document.addEventListener('DOMContentLoaded', function() {
    // "Download txt" 버튼 이벤트 리스너
    document.querySelector('.button.download-txt.valid').addEventListener('click', function() {
        // 로컬 스토리지에서 SMILES 배열 로드
        let smilesArray = JSON.parse(localStorage.getItem('smilesArray'));

        if (!smilesArray || smilesArray.length === 0) {
            alert("No SMILES strings to download.");
            return;
        }

        // SMILES 문자열을 줄바꿈 문자로 결합하여 파일 내용 생성
        let fileContent = smilesArray.join('\n');

        // 텍스트 파일 생성 및 다운로드
        let blob = new Blob([fileContent], {type: 'text/plain'});
        let url = URL.createObjectURL(blob);

        // 다운로드 링크 생성 및 클릭 이벤트 트리거
        let downloadLink = document.createElement('a');
        downloadLink.href = url;
        downloadLink.download = 'valid_smiles.txt'; // 다운로드될 파일명
        document.body.appendChild(downloadLink);
        downloadLink.click();
        document.body.removeChild(downloadLink);
    });

});


document.addEventListener('DOMContentLoaded', function() {
    // .card-box 내부에서 발생하는 클릭 이벤트를 캐치
    document.querySelector('.card-box').addEventListener('click', function(event) {
        // 클릭된 요소가 'Detail' 버튼인지 확인
        if (event.target.classList.contains('detail')) {
            const cardElement = event.target.closest('.card');
            if (cardElement) {
                // .card-contents 내의 첫 번째 h4 태그에서 SMILES 문자열을 추출
                const smiles = cardElement.querySelector('.card-contents h4').textContent;

                // 상세 정보 페이지로 리디렉션합니다.
                window.location.href = `/evaluation-result?smiles=${encodeURIComponent(smiles)}`;
            }
        }
    });
});
