window.onload = function() {
    // Removed MarvinJS initialization
};

document.getElementById('exampleButton').addEventListener('click', function() {
    document.getElementById('smilesInput').value = 'C1=CC=CC=C1O'; // Example SMILES
});

document.getElementById('analyzeButton').addEventListener('click', function() {
    let smiles = document.getElementById('smilesInput').value;
    if (smiles) {
        performADMETPrediction(smiles);
    } else {
        alert("Please enter a SMILES string.");
    }
});


const propertyColors = {
    'Physicochemical Property': '#92a8d1', /* 파스텔 푸른색 */
    'Medicinal Chemistry': '#6b5b95', /* 파스텔 보라색 */
    'Toxicity': '#88b04b', /* 파스텔 녹색 */
    'Absorption': '#f7cac9', /* 파스텔 분홍색 */
    'Distribution': '#ff6f61', /* 파스텔 주황색 */
    'Metabolism': '#779ecb' /* 파스텔 남색 */
};

const propertiesWithFeatures = {
    'Physicochemical Property': [
        'Molecular Weight (MW)',
        'Volume',
        'Density'
    ],
    'Medicinal Chemistry': [
        'QED',
        'SAscore'
    ],
    'Toxicity': [
        'IGC50',
        'LC50FM',
        'LC50DM',
        'Skin Sensitization',
        'Carcinogencity'
    ],
    'Absorption': [
        'Caco-2 Permeability',
        'MDCK Permeability'
    ],
    'Distribution': [
        'PPB',
        'VD'
    ],
    'Metabolism':[
        'CYP1A2 inhibitor',
        'CYP1A2 substrate'
    ]
};


function performADMETPrediction(smiles) {
    // 예시 이미지 경로 (실제 경로로 변경해야 함)
    var imageUrl = "path_to_result_image.png";
    var graphUrl = "path_to_result_graph.png";

    // 체크된 체크박스에서 선택된 속성을 가져옵니다.
    const checkedBoxes = document.querySelectorAll('#propertyFeaturesContainer input[type="checkbox"]:checked');
    const selectedProperties = [];
    checkedBoxes.forEach(cb => {
        if (!cb.id.includes('selectAll')) { // 'Select All' 제외
            const [category, feature] = cb.id.split('-');
            // 임의의 값을 생성합니다.
            const value = (Math.random() * 10).toFixed(3) + ' mg/L';
            selectedProperties.push({ category, feature, value });
        }
    });

    // 결과 이미지와 그래프, 그리고 목업 결과를 표시합니다.
    displayResultImages(imageUrl, graphUrl);
    displayMockResults(selectedProperties);
}

function displayResultImages(imageUrl, graphUrl) {
    const imageContainer = document.querySelector('.image-container');
    const graphContainer = document.querySelector('.graph-container');

    // Assuming you are using placeholders, replace them with actual images
    imageContainer.innerHTML = `<img src="${imageUrl}" alt="Result Image" style="max-width:100%;">`;
    graphContainer.innerHTML = `<img src="${graphUrl}" alt="Result Graph" style="max-width:100%;">`;

    // Remove the 'hidden' class if you have it on these containers to make them visible
    imageContainer.classList.remove('hidden');
    graphContainer.classList.remove('hidden');
}

function createPropertyFeatures() {
    const propertyFeaturesContainer = document.getElementById('propertyFeaturesContainer');
    propertyFeaturesContainer.innerHTML = ''; // Clear previous content

    Object.keys(propertiesWithFeatures).forEach(property => {
        const propertyDiv = document.createElement('div');
        propertyDiv.className = 'property-feature-container';
        propertyDiv.style.borderLeft = `3px solid ${propertyColors[property]}`; // Apply color border

        const header = document.createElement('div');
        header.className = 'property-feature-header';
        header.textContent = property;
        propertyDiv.appendChild(header);

        // 전체 선택 체크박스 생성
        const selectAllCheckbox = document.createElement('input');
        selectAllCheckbox.type = 'checkbox';
        selectAllCheckbox.id = `selectAll-${property}`;
        selectAllCheckbox.className = 'select-all-checkbox';

        const selectAllLabel = document.createElement('label');
        selectAllLabel.htmlFor = `selectAll-${property}`;
        selectAllLabel.textContent = 'Select All';

        const selectAllWrapper = document.createElement('div');
        selectAllWrapper.className = 'select-all-wrapper';
        selectAllWrapper.appendChild(selectAllCheckbox);
        selectAllWrapper.appendChild(selectAllLabel);
        header.appendChild(selectAllWrapper);

        const featureList = document.createElement('ul');
        featureList.className = 'property-feature-list';

        propertiesWithFeatures[property].forEach(feature => {
            const featureItem = document.createElement('li');
            featureItem.className = 'property-feature-item';

            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.id = `${property}-${feature}`;
            checkbox.value = `${property}-${feature}`;

            const label = document.createElement('label');
            label.htmlFor = `${property}-${feature}`;
            label.textContent = feature;

            featureItem.appendChild(checkbox);
            featureItem.appendChild(label);
            featureList.appendChild(featureItem);
        });

        propertyDiv.appendChild(featureList);
        propertyFeaturesContainer.appendChild(propertyDiv);

        // 전체 선택 체크박스 이벤트 리스너 추가
        selectAllCheckbox.addEventListener('change', (e) => {
            const checkboxes = featureList.querySelectorAll('input[type="checkbox"]');
            checkboxes.forEach(checkbox => {
                checkbox.checked = e.target.checked;
            });
        });
    });
}

function displayMockResults(selectedProperties) {
    const resultsContainer = document.querySelector('.results-container');
    resultsContainer.innerHTML = ''; // 이전 결과를 초기화

    // 다운로드 버튼 컨테이너를 생성합니다.
    const buttonsContainer = document.createElement('div');
    buttonsContainer.className = 'buttons-container';

    // CSV 다운로드 버튼을 생성하고 추가합니다.
    const csvButton = document.createElement('button');
    csvButton.textContent = 'Download CSV';
    csvButton.className = 'button-style';
    csvButton.onclick = () => downloadCSV(selectedProperties);
    buttonsContainer.appendChild(csvButton);

    // PDF 다운로드 버튼 생성 및 추가
    // const pdfButton = document.createElement('button');
    // pdfButton.textContent = 'Download PDF';
    // pdfButton.className = 'button-style';
    // pdfButton.onclick = () => downloadPDF(selectedProperties);
    // buttonsContainer.appendChild(pdfButton);

    // 버튼 컨테이너를 결과 컨테이너에 추가합니다.
    resultsContainer.appendChild(buttonsContainer);

    // 카테고리별로 결과를 그룹화하여 결과를 표시합니다.
    selectedProperties.forEach(({ category, feature, value }) => {
        const categoryDiv = resultsContainer.querySelector(`.category-container[data-category="${category}"]`) || createCategoryDiv(category);

        const resultItem = document.createElement('li');
        resultItem.className = 'result-item';
        resultItem.textContent = `${feature}: ${value}`;
        categoryDiv.querySelector('.results-list').appendChild(resultItem);

        if (!categoryDiv.isConnected) {
            resultsContainer.appendChild(categoryDiv);
        }
    });

    resultsContainer.classList.remove('hidden');
}

// 카테고리별 컨테이너를 생성하는 헬퍼 함수
function createCategoryDiv(category) {
    const categoryDiv = document.createElement('div');
    categoryDiv.className = 'category-container';
    categoryDiv.dataset.category = category;
    categoryDiv.style.borderLeft = `3px solid ${propertyColors[category]}`;

    const categoryHeader = document.createElement('div');
    categoryHeader.className = 'category-header';
    categoryHeader.textContent = category;
    categoryDiv.appendChild(categoryHeader);

    const resultsList = document.createElement('ul');
    resultsList.className = 'results-list';
    categoryDiv.appendChild(resultsList);

    return categoryDiv;
}

function downloadCSV(selectedProperties) {
    // CSV 문자열로 변환
    let csvContent = "data:text/csv;charset=utf-8,";
    csvContent += "Category,Property,Value\n"; // CSV 헤더

    selectedProperties.forEach(property => {
        // 카테고리, 프로퍼티, 값 형식으로 CSV에 각 항목 추가
        csvContent += `${property.category},${property.feature},${property.value}\n`;
    });

    // URI 인코딩 후 CSV 내용을 'href'에 할당
    const encodedUri = encodeURI(csvContent);
    // 다운로드 링크 생성
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "ADMET_results.csv");
    document.body.appendChild(link); // 필요하다면 숨겨진 요소로 추가

    link.click(); // 링크를 클릭하여 다운로드 실행
}



function downloadPDF(selectedProperties) {
    // Use the correct syntax to instantiate jsPDF
    const { jsPDF } = window.jspdf.umd;
    const doc = new jsPDF();

    doc.text("Results", 10, 10);
    let y = 20;

    // Loop through selectedProperties, excluding any 'Select All' entries
    selectedProperties.forEach(property => {
        if (property && !property.feature.includes('selectAll')) { // Ensure the property is not undefined and not 'Select All'
            const mockValue = (Math.random() * 10).toFixed(3);
            doc.text(`${property.category} - ${property.feature}: ${mockValue} mg/L`, 10, y);
            y += 10;
        }
    });

    doc.save("results.pdf");
}




window.onload = function() {

    createPropertyFeatures();
};