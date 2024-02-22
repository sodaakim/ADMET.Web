document.addEventListener("DOMContentLoaded", function() {
    // PDF 다운로드 이벤트 리스너는 이미 존재함

    // CSV 다운로드 버튼에 대한 이벤트 리스너 추가
    document.getElementById("download-csv").addEventListener("click", function() {
        downloadCSV(); // CSV 다운로드 함수 실행
    });
});


function convertToCSV(objArray) {
    const array = typeof objArray !== 'object' ? JSON.parse(objArray) : objArray;
    let str = '';

    for (let i = 0; i < array.length; i++) {
        let line = '';
        for (let index in array[i]) {
            if (line !== '') line += ',';

            line += array[i][index];
        }

        str += line + '\r\n';
    }

    return str;
}

function downloadCSV() {
    const data = extractPageData_csv(); // 이전에 정의된 함수를 사용하여 데이터를 추출
    let csvContent = "data:text/csv;charset=utf-8,";

    // 헤더 추가
    csvContent += "Molecule Name,SMILES Notation,Property,Model 1,Model 2,Model 3,Color\n";

    // 데이터 추가
    data.properties.forEach(property => {
        property.contents.forEach(content => {
            csvContent += `${data.molecule_name},${data.smiles_notation},${property.title},${content.model1},${content.model2},${content.model3},${content.color}\n`;
        });
    });

    // CSV 파일 다운로드
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "molecule_data.csv");
    document.body.appendChild(link); // 필요한 경우 Firefox에서 동작하도록 DOM에 링크 추가

    link.click(); // 다운로드 링크 실행
    document.body.removeChild(link); // 다운로드 후 링크 제거
}


function extractPageData_csv() {
    // 페이지에서 데이터 추출 로직 구현
    return {
        molecule_name: document.querySelector('.result-header h2').textContent,
        molecule_image_url: document.querySelector('.molecule-image img').src,
        smiles_notation: document.querySelector('.molecule-smiles h3').textContent,
        properties: Array.from(document.querySelectorAll('.molecule-card .card')).map(card => ({
            title: card.querySelector('.card-title').textContent,
            contents: Array.from(card.querySelectorAll('.card-content')).map(content => ({
                subtitle: content.querySelector('.card-subtitle')?.textContent,
                model1: content.querySelector('.card-num:nth-child(2)')?.textContent,
                model2: content.querySelector('.card-num:nth-child(3)')?.textContent,
                model3: content.querySelector('.card-num:nth-child(4)')?.textContent,
                color: content.querySelector('.color-indicator')?.style.backgroundColor,
            }))
        }))
    };
}
