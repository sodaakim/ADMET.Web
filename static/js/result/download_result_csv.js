
document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('download-result').addEventListener('click', downloadCSV);
});

function downloadCSV() {

    let csvContent = "data:text/csv;charset=utf-8,SMILES";

    selectedProperties.forEach(property => {
        if (property !== 'SMILES') {
            csvContent += `,${property}`;
        }
    });
    csvContent += "\n";

    analysisResults.forEach(result => {
        let row = [result['SMILES']];
        selectedProperties.forEach(property => {
            if (property !== 'SMILES') {
                row.push(result[property] || '');
            }
        });
        csvContent += row.join(",") + "\n";
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


