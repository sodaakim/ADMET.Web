document.addEventListener("DOMContentLoaded", function() {
    document.getElementById("download-pdf").addEventListener("click", function() {
        // 페이지 내용을 JSON으로 변환
        const dataToSend = extractPageData_pdf();

        // 서버로 JSON 데이터 전송
        fetch('/download-pdf', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(dataToSend)
        })
        .then(response => response.blob())
        .then(blob => {
            // PDF 파일로 다운로드
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'molecule_report.pdf';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
        });
    });
});

function extractPageData_pdf() {
    // 페이지에서 데이터 추출 로직 구현
    return {
        molecule_name: document.querySelector('.result-header h2').textContent,
        molecule_image_url: document.querySelector('.molecule-image img').src,
        smiles_notation: document.querySelector('.molecule-smiles h3').textContent,
        properties: Array.from(document.querySelectorAll('.molecule-card .card')).map(card => ({
            title: card.querySelector('.card-title').textContent,
            contents: Array.from(card.querySelectorAll('.card-content.contents')).map(content => ({
                subtitle: content.querySelector('.card-subtitle')?.textContent,
                model1: content.querySelector('.card-num:nth-child(2)')?.textContent,
                model2: content.querySelector('.card-num:nth-child(3)')?.textContent,
                model3: content.querySelector('.card-num:nth-child(4)')?.textContent,
                color: content.querySelector('.color-indicator')?.style.backgroundColor,
            }))
        }))
    };
}
