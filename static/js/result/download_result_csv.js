
document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('download-result').addEventListener('click', downloadCSV);
});

function downloadCSV() {
    let csvContent = "data:text/csv;charset=utf-8,Molecule Name,SMILES";

    const propertiesTitles = Object.keys(analysisResults[0]);
    propertiesTitles.forEach((title) => {
        if (title !== 'SMILES') {
            csvContent += `,${title}`;
        }
    });

    // 데이터 행 추가
    analysisResults.forEach((result) => {
        csvContent += `\n${result['Molecule Name']},${result['SMILES']}`;
        propertiesTitles.forEach((title) => {
            if (title !== 'SMILES') {
                csvContent += `,${result[title]}`;
            }
        });
    });

    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "molecule_properties.csv");
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

