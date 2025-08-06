/* Skyborn Documentation - LaTeX Formula Manager */

// Wait for DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {

    // Track processed formulas to avoid duplicates
    const processedFormulas = new Set();

    // ========== LaTeX Formula Enhancement ==========
    function enhanceLatexFormulas() {
        // Wait for MathJax to render formulas
        setTimeout(() => {
            const formulas = document.querySelectorAll('.MathJax, .math, [class*="MathJax"], mjx-container');

            formulas.forEach(formula => {
                // Skip if already processed or is an image element
                if (processedFormulas.has(formula) ||
                    formula.parentElement.classList.contains('formula-wrapper') ||
                    formula.tagName === 'IMG') {
                    return;
                }

                // Mark as processed
                processedFormulas.add(formula);

                // Create wrapper with enhanced styling
                const wrapper = document.createElement('div');
                wrapper.className = 'formula-wrapper enhanced-formula';

                // Insert wrapper before formula
                formula.parentNode.insertBefore(wrapper, formula);

                // Move formula into wrapper
                wrapper.appendChild(formula);

                // Add click event for zoom
                wrapper.addEventListener('click', (e) => {
                    e.preventDefault();
                    e.stopPropagation();
                    openFormulaModal(formula);
                });

                // Add hover effects
                wrapper.addEventListener('mouseenter', () => {
                    wrapper.style.transform = 'scale(1.05)';
                    wrapper.style.boxShadow = '0 8px 25px rgba(0, 0, 0, 0.15)';
                });

                wrapper.addEventListener('mouseleave', () => {
                    wrapper.style.transform = 'scale(1)';
                    wrapper.style.boxShadow = '0 4px 15px rgba(0, 0, 0, 0.1)';
                });
            });
        }, 1000); // Wait for MathJax to finish rendering
    }

    // ========== Formula Modal with Drag Support ==========
    function openFormulaModal(originalFormula) {
        // Prevent multiple modals
        if (document.querySelector('.formula-modal')) return;

        const modal = createFormulaModal('formula-modal');

        // Clone the formula element properly for display
        const formulaClone = originalFormula.cloneNode(true);

        const content = document.createElement('div');
        content.className = 'modal-formula-content';
        content.appendChild(formulaClone);

        modal.querySelector('.modal-content').appendChild(content);

        // Enhanced styling for formula modal
        content.style.cssText = `
            padding: 40px;
            background: white;
            border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            max-width: 90vw;
            max-height: 80vh;
            overflow: auto;
            text-align: center;
            transform: scale(1.8);
            cursor: move;
            user-select: none;
            position: relative;
        `;

        // Add drag functionality
        makeDraggable(content);

        // Add zoom functionality
        makeZoomable(content);

        document.body.appendChild(modal);

        // Animation
        setTimeout(() => {
            modal.style.opacity = '1';
            content.style.transform = 'scale(1.8) translateY(0)';
        }, 10);
    }

    // ========== Drag Functionality ==========
    function makeDraggable(element) {
        let isDragging = false;
        let currentX = 0;
        let currentY = 0;
        let initialX = 0;
        let initialY = 0;
        let xOffset = 0;
        let yOffset = 0;

        element.addEventListener('mousedown', dragStart);
        document.addEventListener('mousemove', drag);
        document.addEventListener('mouseup', dragEnd);

        // Touch events for mobile
        element.addEventListener('touchstart', dragStart);
        document.addEventListener('touchmove', drag);
        document.addEventListener('touchend', dragEnd);

        function dragStart(e) {
            if (e.type === "touchstart") {
                initialX = e.touches[0].clientX - xOffset;
                initialY = e.touches[0].clientY - yOffset;
            } else {
                initialX = e.clientX - xOffset;
                initialY = e.clientY - yOffset;
            }

            if (e.target === element || element.contains(e.target)) {
                isDragging = true;
                element.style.cursor = 'grabbing';
            }
        }

        function drag(e) {
            if (isDragging) {
                e.preventDefault();

                if (e.type === "touchmove") {
                    currentX = e.touches[0].clientX - initialX;
                    currentY = e.touches[0].clientY - initialY;
                } else {
                    currentX = e.clientX - initialX;
                    currentY = e.clientY - initialY;
                }

                xOffset = currentX;
                yOffset = currentY;

                setTranslate(currentX, currentY, element);
            }
        }

        function setTranslate(xPos, yPos, el) {
            el.style.transform = `scale(1.8) translate(${xPos}px, ${yPos}px)`;
        }

        function dragEnd() {
            if (isDragging) {
                initialX = currentX;
                initialY = currentY;
                isDragging = false;
                element.style.cursor = 'move';
            }
        }
    }

    // ========== Modal Creation Utility ==========
    function createFormulaModal(className) {
        const modal = document.createElement('div');
        modal.className = `modal ${className}`;
        modal.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0, 0, 0, 0.8);
            display: flex;
            justify-content: center;
            align-items: center;
            z-index: 10000;
            opacity: 0;
            transition: opacity 0.3s ease;
            backdrop-filter: blur(5px);
        `;

        const content = document.createElement('div');
        content.className = 'modal-content';
        content.style.cssText = `
            position: relative;
            transform: scale(0.8) translateY(50px);
            transition: all 0.3s cubic-bezier(0.34, 1.56, 0.64, 1);
        `;

        // Close button
        const closeBtn = document.createElement('button');
        closeBtn.innerHTML = '×';
        closeBtn.className = 'modal-close';
        closeBtn.style.cssText = `
            position: absolute;
            top: -20px;
            right: -20px;
            width: 50px;
            height: 50px;
            border: none;
            background: rgba(255, 255, 255, 0.9);
            color: #333;
            font-size: 30px;
            font-weight: bold;
            border-radius: 50%;
            cursor: pointer;
            z-index: 10001;
            display: flex;
            align-items: center;
            justify-content: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            transition: all 0.2s ease;
        `;

        closeBtn.addEventListener('mouseenter', () => {
            closeBtn.style.transform = 'scale(1.1)';
            closeBtn.style.background = '#ff4757';
            closeBtn.style.color = 'white';
        });

        closeBtn.addEventListener('mouseleave', () => {
            closeBtn.style.transform = 'scale(1)';
            closeBtn.style.background = 'rgba(255, 255, 255, 0.9)';
            closeBtn.style.color = '#333';
        });

        // Drag instruction
        const dragHint = document.createElement('div');
        // dragHint.innerHTML = '💡 可以拖拽移动公式';
        dragHint.style.cssText = `
            position: absolute;
            top: -60px;
            left: 50%;
            transform: translateX(-50%);
            background: rgba(255, 255, 255, 0.9);
            color: #333;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 14px;
            font-weight: 500;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.1);
            white-space: nowrap;
            animation: fadeInOut 4s ease-in-out;
        `;

        // Close modal events
        const closeModal = () => {
            modal.style.opacity = '0';
            content.style.transform = 'scale(0.8) translateY(50px)';
            setTimeout(() => modal.remove(), 300);
        };

        closeBtn.addEventListener('click', closeModal);
        modal.addEventListener('click', (e) => {
            if (e.target === modal) closeModal();
        });

        // ESC key to close
        const escHandler = (e) => {
            if (e.key === 'Escape') {
                closeModal();
                document.removeEventListener('keydown', escHandler);
            }
        };
        document.addEventListener('keydown', escHandler);

        content.appendChild(closeBtn);
        content.appendChild(dragHint);
        modal.appendChild(content);

        return modal;
    }

    // ========== CSS Injection ==========
    function injectFormulaStyles() {
        const styles = `
            <style id="latex-formula-styles">
            /* Enhanced Formula Styling */
            .enhanced-formula {
                display: inline-block;
                margin: 15px 10px;
                padding: 20px 25px;
                background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
                border-radius: 15px;
                box-shadow: 0 4px 15px rgba(0, 0, 0, 0.1);
                border: 2px solid rgba(255, 255, 255, 0.8);
                cursor: pointer;
                transition: all 0.3s cubic-bezier(0.34, 1.56, 0.64, 1);
                position: relative;
                overflow: hidden;
            }

            /* 防止MathJax元素被其他脚本处理 */
            .enhanced-formula .MathJax,
            .enhanced-formula mjx-container,
            .enhanced-formula [class*="MathJax"] {
                pointer-events: none !important;
                cursor: inherit !important;
            }

            /* 确保公式包装器内的图片不被通用图片处理器选中 */
            .formula-wrapper img,
            .enhanced-formula img {
                pointer-events: none !important;
                cursor: inherit !important;
            }

            .enhanced-formula::before {
                content: '';
                position: absolute;
                top: 0;
                left: -100%;
                width: 100%;
                height: 100%;
                background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.4), transparent);
                transition: left 0.5s;
            }

            .enhanced-formula:hover::before {
                left: 100%;
            }

            .enhanced-formula:active {
                transform: scale(0.98);
            }

            /* Modal Animations */
            .modal-formula-content {
                animation: modalSlideIn 0.5s cubic-bezier(0.34, 1.56, 0.64, 1);
            }

            @keyframes modalSlideIn {
                from {
                    opacity: 0;
                    transform: scale(0.5) translateY(100px) rotateX(30deg);
                }
                to {
                    opacity: 1;
                    transform: scale(1.8) translateY(0) rotateX(0deg);
                }
            }

            @keyframes fadeInOut {
                0%, 100% { opacity: 0; }
                20%, 80% { opacity: 1; }
            }

            /* Responsive Design */
            @media (max-width: 768px) {
                .enhanced-formula {
                    margin: 10px 5px;
                    padding: 15px;
                }

                .modal-formula-content {
                    padding: 25px !important;
                    transform: scale(1.4) !important;
                }

                @keyframes modalSlideIn {
                    from {
                        opacity: 0;
                        transform: scale(0.5) translateY(100px) rotateX(30deg);
                    }
                    to {
                        opacity: 1;
                        transform: scale(1.4) translateY(0) rotateX(0deg);
                    }
                }
            }

            /* Dark mode support */
            @media (prefers-color-scheme: dark) {
                .enhanced-formula {
                    background: linear-gradient(135deg, #2d3748 0%, #1a202c 100%);
                    border-color: rgba(255, 255, 255, 0.1);
                    box-shadow: 0 4px 15px rgba(0, 0, 0, 0.3);
                }

                .modal-close {
                    background: rgba(45, 55, 72, 0.9) !important;
                    color: white !important;
                }
            }
            </style>
        `;

        // Remove existing styles to avoid conflicts
        const existingStyles = document.getElementById('latex-formula-styles');
        if (existingStyles) {
            existingStyles.remove();
        }

        document.head.insertAdjacentHTML('beforeend', styles);
    }

    // ========== Initialize Formula Enhancements ==========
    function initializeFormulaEnhancements() {
        injectFormulaStyles();
        enhanceLatexFormulas();

        // Re-enhance when new content is loaded (for dynamic content)
        const observer = new MutationObserver((mutations) => {
            mutations.forEach((mutation) => {
                if (mutation.addedNodes.length > 0) {
                    setTimeout(() => {
                        enhanceLatexFormulas();
                    }, 500);
                }
            });
        });

        observer.observe(document.body, {
            childList: true,
            subtree: true
        });
    }

    // Start initialization
    initializeFormulaEnhancements();

    // Re-initialize after MathJax finishes (if present)
    if (window.MathJax) {
        window.MathJax.startup.promise.then(() => {
            setTimeout(() => {
                processedFormulas.clear(); // Clear processed formulas
                enhanceLatexFormulas();
            }, 500);
        });
    }

    console.log('✅ LaTeX Formula Manager loaded successfully!');
});
