document.addEventListener("DOMContentLoaded", function() {
    document.getElementById("download-csv").addEventListener("click", function() {
        downloadCSV();
    });
});

function downloadCSV() {
    const data = extractPageData_csv();
    let csvContent = "data:text/csv;charset=utf-8,";

    csvContent += "Molecule Name,SMILES Notation,Property,Property,FPADMET,OPERA,SSBIO\n";

    data.properties.forEach((property, index) => {
        property.contents.forEach((content, contentIndex) => {
            let moleculeRow = index === 0 && contentIndex === 0 ? `${data.molecule_name},${data.smiles_notation},` : ',,';
            let titleRow = contentIndex === 0 ? `${property.title},` : ',';
            csvContent += `${moleculeRow}${titleRow}${content.subtitle},${content.model1},${content.model2},${content.model3}\n`;
        });
    });


    // CSV 다운로드
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
        molecule_name: document.querySelector('.result-header h2').textContent,
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
