window.onload = adjustContainerHeight;
window.onresize = adjustContainerHeight;

function adjustContainerHeight() {
    // 카드들의 부모 컨테이너인 `.molecule-card` 선택
    const moleculeCard = document.querySelector('.molecule-card');
    // `.molecule-card` 내의 모든 `.card` 요소의 높이 합계 계산
    const totalHeight = Array.from(moleculeCard.children).reduce((total, card) => {
        return total + card.offsetHeight;
    }, 0);


    // `.molecule-title` 요소의 높이 계산
    const moleculeTitle = document.querySelector('.molecule-title');
    const titleHeight = moleculeTitle.offsetHeight;

    // 계산된 전체 높이의 절반 값
    const halfHeight = totalHeight / 2;

    moleculeCard.style.maxHeight = `${halfHeight+200}px`;

    // 계산된 전체 카드 높이의 절반 값에 `.molecule-title`의 높이를 더함
    const adjustedHeight = (totalCardHeight / 2) + titleHeight;

    const resultBox = document.querySelector('.result-box');
    resultBox.style.height = `${adjustedHeight}px`;
    const moleculeContents = document.querySelector('.molecule-contents');
    moleculeContents.style.height = `${halfHeight}px`;
}


document.querySelector('.tooltip').addEventListener('mouseover', function(event) {
  var tooltipText = this.querySelector('.tooltiptext');
  tooltipText.style.visibility = 'visible';
});

document.querySelector('.tooltip').addEventListener('mouseout', function(event) {
  var tooltipText = this.querySelector('.tooltiptext');
  tooltipText.style.visibility = 'hidden';
});
