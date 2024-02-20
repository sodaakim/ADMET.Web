document.addEventListener('DOMContentLoaded', function() {
  document.querySelectorAll('.toggle-button').forEach(button => {
    button.addEventListener('click', function() {
      const card = this.closest('.property-feature-card');
      const cardContent = card.querySelector('.card-content');
      const contentSeparator = card.querySelector('.content-separator');
      if (cardContent.style.display === 'none') {
        cardContent.style.display = 'block';
        contentSeparator.style.display = 'block';
        this.textContent = '-'; // 토글 상태 아이콘으로 대체 가능
      } else {
        cardContent.style.display = 'none';
        contentSeparator.style.display = 'none';
        this.textContent = '+'; // 토글 상태 아이콘으로 대체 가능
      }
    });
  });

  // Select All 기능 구현
  document.querySelectorAll('.select-all-checkbox').forEach(checkbox => {
    checkbox.addEventListener('change', function() {
      const card = this.closest('.property-feature-card');
      const checkboxes = card.querySelectorAll('.property-feature-item input[type="checkbox"]');
      checkboxes.forEach(cb => cb.checked = this.checked);
    });
  });
});
