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

document.addEventListener('DOMContentLoaded', function() {
    var toggleButton = document.querySelector('.toggleDownloadOptions');
    var container = document.getElementById('propertyFeaturesContainer');

    // Ensure the container is hidden and the button text is set to '+' when the page loads
    container.style.display = 'none'; // Make sure the container is hidden
    toggleButton.textContent = '+'; // Set the button text to '+'

    // Add click event listener to the button
    toggleButton.addEventListener('click', function() {
        // Toggle the display style
        var isHidden = container.style.display === 'none';
        container.style.display = isHidden ? 'grid' : 'none';
        // Toggle the button content
        toggleButton.textContent = isHidden ? '−' : '+';
    });
});

