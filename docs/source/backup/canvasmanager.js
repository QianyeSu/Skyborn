/**
 * Canvas manager for handling canvas setup, resizing, and context management.
 * @module canvasManager
 */

/**
 * Validates a dimension input.
 * @param {string|number} dim - Dimension value or 'full'.
 * @param {string} name - Parameter name for error messages.
 * @returns {number} Validated dimension.
 * @throws {TypeError} If the dimension is invalid.
 */
function validate(dim, name) {
    if (dim === 'full') {
        return name === 'width' ? window.innerWidth : window.innerHeight;
    }
    const parsed = parseInt(dim, 10);
    if (!isNaN(parsed) && parsed > 0) {
        return parsed;
    }
    throw new TypeError(`${name} must be a positive number or 'full'`);
}

/**
 * Attaches canvas manager to a canvas element.
 * @param {HTMLCanvasElement|string} elementOrId - Canvas element or ID.
 * @param {Object} [config] - Optional configuration object.
 * @param {string|number} [config.width='full'] - Initial width or 'full'.
 * @param {string|number} [config.height='full'] - Initial height or 'full'.
 * @param {'2d'|'webgl'|'webgl2'} [config.context='2d'] - Context type.
 * @param {Object} [config.contextOptions={}] - Context creation options.
 * @param {boolean} [config.autoInit=false] - Whether to automatically initialize context and size.
 * @param {boolean} [config.autoListen=false] - Whether to automatically listen for resize events.
 * @param {number} [config.resizeDebounce=250] - Debounce time for resize events (ms).
 * @param {boolean} [config.listenDpr=false] - Whether to listen for device pixel ratio changes.
 * @param {number} [config.dprInterval=250] - Interval for DPR checking (ms).
 * @returns {CanvasController} Canvas controller instance.
 * @throws {TypeError} If the argument is invalid.
 */
function attach(elementOrId, config = {}) {
    if (!elementOrId) {
        throw new TypeError('attach requires a canvas element or string ID');
    }

    const el = typeof elementOrId === 'string'
        ? document.getElementById(elementOrId)
        : elementOrId;

    if (!(el instanceof HTMLCanvasElement)) {
        throw new TypeError('attach expects a canvas element or string ID');
    }

    // Extract config options with defaults
    const {
        width: configWidth = 'full',
        height: configHeight = 'full',
        context: contextType = '2d',
        contextOptions = {},
        autoInit = false,
        autoListen = false,
        resizeDebounce = 250,
        listenDpr = false,
        dprInterval = 250
    } = config;

    /** @type {CanvasRenderingContext2D|WebGLRenderingContext|WebGL2RenderingContext|null} */
    let ctx = null;
    let width = 0;
    let height = 0;
    let dpr = window.devicePixelRatio || 1;
    let resizeTimeout = null;
    let resizeHandler = null;
    let dprIntervalID = null;
    let isDestroyed = false;

    /**
     * Throws if the controller is destroyed.
     * @private
     */
    function checkDestroyed() {
        if (isDestroyed) {
            throw new Error('Canvas manager has been destroyed');
        }
    }

    /**
     * Canvas controller interface.
     * @typedef {Object} CanvasController
     * @property {function(string|number=, string|number=): CanvasController} resize - Resize the canvas.
     * @property {function(string=, Object=): CanvasRenderingContext2D | WebGLRenderingContext | WebGL2RenderingContext} context - Get the canvas context.
     * @property {function(string=, number=): CanvasController} listen - Listen for resize or dpr changes.
     * @property {function(): void} destroy - Destroy the controller and clean up.
     * @property {function(string|number, string): number} validate - Dimension validator.
     * @property {HTMLCanvasElement} el - The canvas element.
     * @property {CanvasRenderingContext2D|WebGLRenderingContext|WebGL2RenderingContext|null} ctx - The canvas context.
     * @property {number} width - The logical width.
     * @property {number} height - The logical height.
     * @property {number} dpr - The device pixel ratio.
     * @property {boolean} isDestroyed - Whether the controller is destroyed.
     */

    const controller = {
        /**
         * Apply DPR transform and optionally clear the canvas.
         * Handles 2D vs WebGL contexts appropriately.
         * @private
         * @param {boolean} clear
         */
        applyTransformAndClear(clear) {
            if (!ctx) return;

            if (ctx instanceof CanvasRenderingContext2D) {
                ctx.setTransform(1, 0, 0, 1, 0, 0);
                if (clear) {
                    ctx.clearRect(0, 0, el.width, el.height);
                }
                ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            }

            else if (ctx instanceof WebGLRenderingContext || ctx instanceof WebGL2RenderingContext) {
                if (clear) {
                    const gl = ctx;
                    gl.viewport(0, 0, el.width, el.height);
                    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
                }
            }
        },

        /**
         * Resize the canvas.
         * @param {string|number} [w='full'] - Width or 'full'.
         * @param {string|number} [h='full'] - Height or 'full'.
         * @param {boolean} [clear=true] - Whether to clear the canvas.
         * @returns {CanvasController}
         */
        resize(w = 'full', h = 'full', clear = true) {
            checkDestroyed();
            width = validate(w, 'width');
            height = validate(h, 'height');
            dpr = window.devicePixelRatio || 1;

            el.width = width * dpr;
            el.height = height * dpr;
            el.style.width = `${width}px`;
            el.style.height = `${height}px`;

            controller.applyTransformAndClear(clear);

            return controller;
        },

        /**
         * Get the canvas context.
         * @param {'2d'|'webgl'|'webgl2'} [type='2d'] - Context type.
         * @param {Object} [options={}] - Context options.
         * @returns {CanvasRenderingContext2D|WebGLRenderingContext|WebGL2RenderingContext}
         * @throws {Error} If context cannot be obtained.
         */
        context(type = '2d', options = {}) {
            checkDestroyed();
            if (!ctx) {
                ctx = el.getContext(type, options);
                if (!ctx) {
                    throw new Error(`Failed to get '${type}' context`);
                }
                if (width === 0 || height === 0) {
                    const w = el.clientWidth;
                    const h = el.clientHeight;
                    if (!w || !h) {
                        console.warn('Canvas has no visible size; falling back to default 300x150.');
                    }
                    controller.resize(w || 300, h || 150);
                }
            }
            return ctx;
        },

        /**
         * Listen for resize or dpr changes.
         * @param {string} [signal='resize'] - 'resize' or 'dpr'.
         * @param {number} [time=250] - Debounce time or interval (ms).
         * @returns {CanvasController}
         * @throws {TypeError} If time is invalid.
         */
        listen(signal = 'resize', time = 250) {
            checkDestroyed();
            if (typeof time !== 'number' || time < 0) {
                throw new TypeError('time must be a non-negative number');
            }

            if (signal === 'resize') {
                if (resizeHandler) {
                    window.removeEventListener('resize', resizeHandler);
                }
                resizeHandler = () => {
                    clearTimeout(resizeTimeout);
                    resizeTimeout = setTimeout(() => {
                        if (!isDestroyed) {
                            controller.resize(width === 0 ? 'full' : width, height === 0 ? 'full' : height);
                        }
                    }, time);
                };
                window.addEventListener('resize', resizeHandler);
            }

            if (signal === 'dpr') {
                if (dprIntervalID) {
                    clearInterval(dprIntervalID);
                }
                dprIntervalID = setInterval(() => {
                    if (isDestroyed) {
                        clearInterval(dprIntervalID);
                        return;
                    }
                    const current = window.devicePixelRatio || 1;
                    if (current !== dpr) {
                        dpr = current;
                        controller.resize(width === 0 ? 'full' : width, height === 0 ? 'full' : height);
                    }
                }, time);
            }

            return controller;
        },

        /**
         * Initialize the canvas with specified context and dimensions.
         * @param {string|number} [w='full'] - Width or 'full'.
         * @param {string|number} [h='full'] - Height or 'full'.
         * @param {'2d'|'webgl'|'webgl2'} [type='2d'] - Context type.
         * @param {Object} [options={}] - Context options.
         * @returns {CanvasController} - The canvas controller instance.
         * @throws {TypeError} If dimensions are invalid or context type is unsupported.
         * @throws {Error} If context cannot be obtained.
         */
        initWithContext(w = 'full', h = 'full', type = '2d', options = {}) {
            checkDestroyed();
            if (type !== '2d' && type !== 'webgl' && type !== 'webgl2') {
                throw new TypeError(`Unsupported context type: ${type}`);
            }
            width = validate(w, 'width');
            height = validate(h, 'height');
            dpr = window.devicePixelRatio || 1;
            el.width = width * dpr;
            el.height = height * dpr;
            el.style.width = `${width}px`;
            el.style.height = `${height}px`;
            ctx = el.getContext(type, options);
            if (!ctx) {
                throw new Error(`Failed to get '${type}' context`);
            }
            controller.applyTransformAndClear(true);
            return controller;
        },

        /**
         * Destroy the controller and clean up.
         */
        destroy() {
            isDestroyed = true;
            if (resizeTimeout) {
                clearTimeout(resizeTimeout);
                resizeTimeout = null;
            }
            if (resizeHandler) {
                window.removeEventListener('resize', resizeHandler);
                resizeHandler = null;
            }
            if (dprIntervalID) {
                clearInterval(dprIntervalID);
                dprIntervalID = null;
            }
            ctx = null;
        },

        validate,
        get el() { checkDestroyed(); return el; },
        get ctx() { checkDestroyed(); return ctx; },
        get width() { return width; },
        get height() { return height; },
        get dpr() { return dpr; },
        get isDestroyed() { return isDestroyed; },
    };

    if (autoInit) {
        controller.initWithContext(configWidth, configHeight, contextType, contextOptions);
    }

    if (autoListen) {
        controller.listen('resize', resizeDebounce);
    }

    if (listenDpr) {
        controller.listen('dpr', dprInterval);
    }

    return controller;
}

/**
 * Canvas manager factory.
 * @returns {{attach: typeof attach, validate: typeof validate}}
 */
function canvasManager() {
    return { attach, validate };
}

export default canvasManager;
