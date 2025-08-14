/**
 * FPSCounter displays frames per second and memory usage in a fixed overlay.
 * @example
 * const fps = new FPSCounter(); // Adds overlay to document.body
 * const fps = new FPSCounter(document.getElementById('container'));
 */
export default class FPSCounter {
    /**
     * @typedef {Object} FPSCounterOptions
     * @property {boolean} [manual=false] - If true, disables internal animation loop; use update() manually.
     */
    /**
     * Creates an FPSCounter instance and attaches it to the given container.
     * @param {HTMLElement} [container=document.body] - The DOM element to attach the counter to.
     * @param {FPSCounterOptions} [options] - Optional settings.
     * @throws {Error} If container is not a valid DOM element.
     */
    constructor(container = document.body, options = {}) {
        if (!(container instanceof HTMLElement)) {
            throw new Error('Container must be a valid DOM element.');
        }

        this.fps = 0;
        this.frameCount = 0;
        this.lastTime = performance.now();
        this.fpsInterval = 1000;
        this.lastFpsUpdate = this.lastTime;
        this.framesSinceLastUpdate = 0;

        this.container = container;
        this.manual = !!options.manual;
        this.bare = !!options.bare;

        if (this.bare) {
            // In bare mode, use the container directly for text updates
            this.fpsElement = container;
            this.frameCountElement = null;
            this.memElement = null;
        } else {
            const wrapper = document.createElement('div');
            wrapper.className = 'fps-counter';
            wrapper.style.cssText = `    position: fixed;
                top: 0;
                left: 0;
                background: rgba(0,0,0,0.7);
                color: lime;
                font: 12px monospace;
                padding: 5px;
                z-index: 9999;`;

            this.frameCountElement = document.createElement('div');
            this.frameCountElement.className = 'fps-frame-count';
            wrapper.appendChild(this.frameCountElement);

            this.fpsElement = document.createElement('div');
            this.fpsElement.className = 'fps-value';
            wrapper.appendChild(this.fpsElement);

            this.memElement = document.createElement('div');
            this.memElement.className = 'fps-memory-usage';
            this.memElement.style.cssText = `    cursor: help;`;
            this.memElement.title = 'Only works in Chrome-ish browsers';
            wrapper.appendChild(this.memElement);

            this.container.appendChild(wrapper);
        }

        if (!this.manual) {
            this.animate();
        }
    }

    /**
     * Internal animation loop that updates FPS and memory usage.
     * @private
     */
    animate() {
        const currentTime = performance.now();
        this.frameCount++;
        this.framesSinceLastUpdate++;

        if (currentTime - this.lastFpsUpdate >= this.fpsInterval) {
            this.fps = Math.round((this.framesSinceLastUpdate * 1000) / (currentTime - this.lastFpsUpdate));
            this.lastFpsUpdate = currentTime;
            this.framesSinceLastUpdate = 0;
            if (this.bare) {
                this.fpsElement.textContent = `${this.fps} FPS`;
            } else {
                this.fpsElement.textContent = `${this.fps} FPS`;
                if (performance.memory && this.memElement) {
                    const { usedJSHeapSize, totalJSHeapSize } = performance.memory;
                    const usedMB = (usedJSHeapSize / 1024 / 1024).toFixed(2);
                    const totalMB = (totalJSHeapSize / 1024 / 1024).toFixed(2);
                    this.memElement.textContent = `Heap: ${usedMB} / ${totalMB} MB`;
                }
            }
        }
        if (!this.bare && this.frameCountElement) {
            this.frameCountElement.textContent = `Frames: ${this.frameCount}`;
        }
        this.lastTime = currentTime;
        if (!this.manual) {
            requestAnimationFrame(() => this.animate());
        }
    }

    /**
     * Manually update the FPS and frame count display. Useful for external simulation loops.
     * @param {number} fps - The current FPS value to display.
     * @param {number} frameCount - The current frame count to display.
     */
    update(fps, frameCount) {
        if (typeof fps === 'number') {
            this.fps = fps;
            this.fpsElement.textContent = `${fps} FPS`;
        }
        if (!this.bare && typeof frameCount === 'number' && this.frameCountElement) {
            this.frameCount = frameCount;
            this.frameCountElement.textContent = `Frames: ${frameCount}`;
        }
        if (!this.bare && performance.memory && this.memElement) {
            const { usedJSHeapSize, totalJSHeapSize } = performance.memory;
            const usedMB = (usedJSHeapSize / 1024 / 1024).toFixed(2);
            const totalMB = (totalJSHeapSize / 1024 / 1024).toFixed(2);
            this.memElement.textContent = `Heap: ${usedMB} / ${totalMB} MB`;
        }
    }

    /**
     * Alias for update() since people seem to think it ought to exist
     * @param {number} fps - The current FPS value to display.
     * @param {number} frameCount - The current frame count to display.
     */
    tick(fps, framecount) {
        this.update(fps, framecount);
    }

}
