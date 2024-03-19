document.addEventListener("DOMContentLoaded", function() {
    document.getElementById("download-csv").addEventListener("click", function() {
        downloadCSV();
    });
});

function downloadCSV() {
    const data = extractPageData_csv();
    let csvContent = "data:text/csv;charset=utf-8,SMILES Notation";

    // 헤더 구성
    data.properties.forEach((property, index) => {
        property.contents.forEach(content => {
            if (index === 0 || index === 1) { // 프로퍼티 1, 2인 경우 Model3만 출력
                csvContent += `,${property.title} - ${content.subtitle} Model3`;
            } else { // 그 외의 경우 Model1, Model2, Model3 모두 출력
                csvContent += `,${property.title} - ${content.subtitle} Model1,${property.title} - ${content.subtitle} Model2,${property.title} - ${content.subtitle} Model3`;
            }
        });
    });
    csvContent += "\n";

    // 데이터 행 구성
    let rowData = data.smiles_notation;
    data.properties.forEach((property, index) => {
        property.contents.forEach(content => {
            if (index === 0 || index === 1) {
                rowData += `,${content.model3}`;
            } else {
                rowData += `,${content.model1},${content.model2},${content.model3}`;
            }
        });
    });
    csvContent += rowData;

    // CSV 다운로드 로직
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "molecule_data.csv");
    document.body.appendChild(link);

    link.click();
    document.body.removeChild(link);
}



function extractPageData_csv() {
    return {
        // molecule_name: document.querySelector('.result-header h2').textContent,
        smiles_notation: document.querySelector('.molecule-smiles h3').textContent,
        properties: Array.from(document.querySelectorAll('.molecule-card .card')).map(card => ({
            title: card.querySelector('.card-title').textContent,
            contents: Array.from(card.querySelectorAll('.card-content.contents')).map(content => ({
                subtitle: content.querySelector('.card-subtitle')?.textContent,
                model1: content.querySelector('.card-num:nth-child(2)')?.textContent,
                model2: content.querySelector('.card-num:nth-child(3)')?.textContent,
                model3: content.querySelector('.card-num:nth-child(4)')?.textContent,
            }))
        }))
    };
}
