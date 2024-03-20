document.addEventListener('DOMContentLoaded', function() {
    document.querySelectorAll('.property-feature-item input[type="checkbox"]').forEach(checkbox => {
        checkbox.checked = true;
    });
    document.querySelectorAll('.select-all-checkbox').forEach(selectAllCheckbox => {
        selectAllCheckbox.checked = true;
    });

    document.querySelectorAll('.toggle-button').forEach(button => {
        button.addEventListener('click', function() {
        const card = this.closest('.property-feature-card');
        const cardContent = card.querySelector('.card-content');
        const contentSeparator = card.querySelector('.content-separator');
        if (cardContent.style.display === 'none') {
            cardContent.style.display = 'block';
            contentSeparator.style.display = 'block';
            this.textContent = '-';
        } else {
            cardContent.style.display = 'none';
            contentSeparator.style.display = 'none';
            this.textContent = '+';
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

document.addEventListener('DOMContentLoaded', function() {
    var toggleButton = document.querySelector('.toggleDownloadOptions');
    var container = document.getElementById('propertyFeaturesContainer');

    container.style.display = 'none';
    toggleButton.textContent = '+';

    toggleButton.addEventListener('click', function() {
        var isHidden = container.style.display === 'none';
        container.style.display = isHidden ? 'grid' : 'none';
        toggleButton.textContent = isHidden ? '−' : '+';
    });
});
