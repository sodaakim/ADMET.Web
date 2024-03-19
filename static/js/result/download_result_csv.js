
document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('download-result').addEventListener('click', downloadCSV);
});

function downloadCSV() {
    let csvContent = "data:text/csv;charset=utf-8,SMILES";

    const propertiesTitles = Object.keys(analysisResults[0]).filter(title => title !== "Molecule Name" && title !== "Image URL");
    propertiesTitles.forEach((title) => {
        if (title !== 'SMILES') {
            csvContent += `,${title}`;
        }
    });

    // 각 결과에 대한 데이터 행 추가
    analysisResults.forEach((result) => {
        csvContent += `\n${result['SMILES']}`;
        propertiesTitles.forEach((title) => {
            if (title !== 'SMILES') {
                csvContent += `,${result[title]}`;
            }
        });
    });

    // CSV 파일 다운로드
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "molecule_properties.csv");
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}


