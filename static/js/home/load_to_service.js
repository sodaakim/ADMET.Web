document.addEventListener('DOMContentLoaded', () => {
    const manualBox = document.querySelector('.manual');

    manualBox.addEventListener('click', () => {
        window.location.href = '/service_manual';
    });
});

document.addEventListener('DOMContentLoaded', () => {
    const manualBox = document.querySelector('.screening');

    manualBox.addEventListener('click', () => {
        window.location.href = '/analyze';
    });
});
