document.addEventListener('DOMContentLoaded', () => {
    function startCounter(className, start, end, intervalTime) {
        let counter = document.querySelector(`.${className}`);
        let count = start;

        counter.textContent = count;

        const interval = setInterval(() => {
            if (count < end) {
                count++;
                counter.textContent = count;
            } else {
                clearInterval(interval);
            }
        }, intervalTime);
    }

    startCounter('counting1', 80, 100, 30);
    startCounter('counting2', 100, 200, 15);
    startCounter('counting3', 200, 300, 10);
});
