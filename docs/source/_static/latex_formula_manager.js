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

                // Check if formula is inline (within a paragraph or text block)
                const isInline = isFormulaInline(formula);

                // Only add wrapper and styling for standalone formulas (not inline ones)
                if (!isInline) {
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
                }
            });
        }, 1000); // Wait for MathJax to finish rendering
    }

    // Function to determine if a formula is inline (within text)
    function isFormulaInline(formula) {
        const parent = formula.parentElement;

        // Check if parent is a paragraph, span, or other text container
        const textContainers = ['p', 'span', 'div', 'td', 'th', 'li', 'dd', 'dt'];
        const parentTag = parent.tagName.toLowerCase();

        // If parent is a text container, check if there's text content around the formula
        if (textContainers.includes(parentTag)) {
            const parentText = parent.textContent.trim();
            const formulaText = formula.textContent.trim();

            // If the parent contains more text than just the formula, it's likely inline
            if (parentText.length > formulaText.length + 10) {
                return true;
            }

            // Check if there are text nodes before or after the formula
            const siblings = Array.from(parent.childNodes);
            const formulaIndex = siblings.indexOf(formula);

            // Check for text content before the formula
            for (let i = 0; i < formulaIndex; i++) {
                const node = siblings[i];
                if (node.nodeType === Node.TEXT_NODE && node.textContent.trim().length > 0) {
                    return true;
                }
                if (node.nodeType === Node.ELEMENT_NODE && node.textContent.trim().length > 0) {
                    return true;
                }
            }

            // Check for text content after the formula
            for (let i = formulaIndex + 1; i < siblings.length; i++) {
                const node = siblings[i];
                if (node.nodeType === Node.TEXT_NODE && node.textContent.trim().length > 0) {
                    return true;
                }
                if (node.nodeType === Node.ELEMENT_NODE && node.textContent.trim().length > 0) {
                    return true;
                }
            }
        }

        // Check if the formula is within a math block (which should be standalone)
        if (parent.classList.contains('math') || parent.classList.contains('displaymath')) {
            return false;
        }

        // Default to standalone (not inline) for safety
        return false;
    }

    // ========== Formula Modal with Image-like Effect ==========
    function openFormulaModal(originalFormula) {
        // Prevent multiple modals
        if (document.querySelector('.formula-modal')) return;

        // Create overlay (similar to image modal)
        const overlay = document.createElement('div');
        overlay.className = 'formula-modal-overlay';
        overlay.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0, 0, 0, 0.9);
            z-index: 10000;
            opacity: 0;
            transition: opacity 0.3s ease;
            display: flex;
            align-items: center;
            justify-content: center;
            backdrop-filter: blur(10px);
        `;

        // Create modal container (similar to image modal)
        const modal = document.createElement('div');
        modal.className = 'formula-modal';
        modal.style.cssText = `
            background: white;
            border-radius: 20px;
            padding: 2rem;
            max-width: 90vw;
            max-height: 90vh;
            min-width: 300px;
            min-height: 200px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.5);
            transform: scale(0.3) rotate(5deg);
            transition: all 0.4s cubic-bezier(0.175, 0.885, 0.32, 1.275);
            opacity: 0;
            position: relative;
            overflow: hidden;
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
        `;

        // Create close button (exactly like image modal)
        const closeBtn = document.createElement('button');
        closeBtn.innerHTML = '×';
        closeBtn.style.cssText = `
            position: absolute;
            top: 10px;
            right: 15px;
            background: rgba(255, 255, 255, 0.9);
            border: none;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            font-size: 1.5rem;
            cursor: pointer;
            color: #333;
            transition: all 0.3s ease;
            z-index: 10001;
            display: flex;
            align-items: center;
            justify-content: center;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.2);
        `;

        closeBtn.addEventListener('mouseenter', () => {
            closeBtn.style.background = 'rgba(255, 0, 0, 0.1)';
            closeBtn.style.color = '#ff0000';
            closeBtn.style.transform = 'scale(1.1)';
        });

        closeBtn.addEventListener('mouseleave', () => {
            closeBtn.style.background = 'rgba(255, 255, 255, 0.9)';
            closeBtn.style.color = '#333';
            closeBtn.style.transform = 'scale(1)';
        });

        // Clone the formula element with enhanced styling (like image)
        const formulaClone = originalFormula.cloneNode(true);
        formulaClone.style.cssText = `
            max-width: 100%;
            max-height: 100%;
            border: none !important;
            border-radius: 10px !important;
            padding: 1rem !important;
            background: none !important;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.3) !important;
            margin: 0 auto;
            display: block;
            animation: formulaZoom 0.5s ease-out;
            font-size: 2rem !important;
            cursor: grab;
            user-select: none;
        `;

        // Add drag and zoom functionality (exactly like image modal)
        let isDragging = false;
        let startX, startY, translateX = 0, translateY = 0;
        let scale = 1;
        let lastTranslateX = 0, lastTranslateY = 0;

        // Mouse down - start dragging
        formulaClone.addEventListener('mousedown', (e) => {
            e.preventDefault();
            isDragging = true;
            startX = e.clientX - translateX;
            startY = e.clientY - translateY;
            formulaClone.style.cursor = 'grabbing';
            formulaClone.style.userSelect = 'none';
            formulaClone.style.transition = 'none';
        });

        // Mouse move - drag formula
        document.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            e.preventDefault();
            translateX = e.clientX - startX;
            translateY = e.clientY - startY;
            formulaClone.style.transform = `translate(${translateX}px, ${translateY}px) scale(${scale})`;
        });

        // Mouse up - end dragging
        document.addEventListener('mouseup', (e) => {
            if (isDragging) {
                isDragging = false;
                formulaClone.style.cursor = 'grab';
                formulaClone.style.userSelect = 'auto';
                formulaClone.style.transition = 'transform 0.2s ease';
                lastTranslateX = translateX;
                lastTranslateY = translateY;
            }
        });

        // Scroll wheel zoom (optimized zoom center point)
        formulaClone.addEventListener('wheel', (e) => {
            e.preventDefault();
            const rect = formulaClone.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            const delta = e.deltaY > 0 ? 0.9 : 1.1;
            const newScale = Math.max(0.5, Math.min(5, scale * delta));

            if (newScale !== scale) {
                const scaleChange = newScale / scale;
                translateX = translateX * scaleChange + mouseX * (1 - scaleChange);
                translateY = translateY * scaleChange + mouseY * (1 - scaleChange);

                scale = newScale;
                formulaClone.style.transform = `translate(${translateX}px, ${translateY}px) scale(${scale})`;

                lastTranslateX = translateX;
                lastTranslateY = translateY;
            }
        });

        // Double click to reset formula position and zoom
        formulaClone.addEventListener('dblclick', (e) => {
            e.preventDefault();
            translateX = 0;
            translateY = 0;
            scale = 1;
            formulaClone.style.transition = 'transform 0.3s ease';
            formulaClone.style.transform = `translate(0px, 0px) scale(1)`;
            setTimeout(() => {
                formulaClone.style.transition = 'transform 0.2s ease';
            }, 300);
            lastTranslateX = 0;
            lastTranslateY = 0;
        });

        // Prevent default drag behavior
        formulaClone.addEventListener('dragstart', (e) => {
            e.preventDefault();
        });

        // Assemble modal content
        modal.appendChild(closeBtn);
        modal.appendChild(formulaClone);
        overlay.appendChild(modal);
        document.body.appendChild(overlay);

        // Prevent page scrolling
        document.body.style.overflow = 'hidden';

        // Animated display (same as image modal)
        requestAnimationFrame(() => {
            overlay.style.opacity = '1';
            modal.style.transform = 'scale(1) rotate(0deg)';
            modal.style.opacity = '1';
        });

        // Close function
        function closeModal() {
            overlay.style.opacity = '0';
            modal.style.transform = 'scale(0.3) rotate(-5deg)';
            modal.style.opacity = '0';

            setTimeout(() => {
                document.body.removeChild(overlay);
                document.body.style.overflow = '';
            }, 400);
        }

        // Bind close events
        closeBtn.addEventListener('click', closeModal);
        overlay.addEventListener('click', (e) => {
            if (e.target === overlay) {
                closeModal();
            }
        });

        // ESC key to close
        const escHandler = (e) => {
            if (e.key === 'Escape') {
                closeModal();
                document.removeEventListener('keydown', escHandler);
            }
        };
        document.addEventListener('keydown', escHandler);

        // Add background effects (similar to image modal)
        addFormulaModalEffects(modal);
    }

    // Add background effects for formula modal (similar to image modal)
    function addFormulaModalEffects(modal) {
        const canvas = document.createElement('canvas');
        canvas.style.cssText = `
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            pointer-events: none;
            z-index: -1;
            opacity: 0.3;
        `;

        modal.appendChild(canvas);

        const ctx = canvas.getContext('2d');
        const dots = [];

        function resizeCanvas() {
            canvas.width = modal.offsetWidth;
            canvas.height = modal.offsetHeight;
        }
        resizeCanvas();

        // Create decorative dot array
        for (let i = 0; i < 20; i++) {
            dots.push({
                x: Math.random() * canvas.width,
                y: Math.random() * canvas.height,
                size: Math.random() * 3 + 1,
                speedX: (Math.random() - 0.5) * 0.5,
                speedY: (Math.random() - 0.5) * 0.5,
                opacity: Math.random() * 0.5 + 0.2
            });
        }

        function animate() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            dots.forEach(dot => {
                dot.x += dot.speedX;
                dot.y += dot.speedY;

                if (dot.x < 0 || dot.x > canvas.width) dot.speedX *= -1;
                if (dot.y < 0 || dot.y > canvas.height) dot.speedY *= -1;

                ctx.beginPath();
                ctx.arc(dot.x, dot.y, dot.size, 0, Math.PI * 2);
                ctx.fillStyle = `rgba(100, 150, 255, ${dot.opacity})`;
                ctx.fill();
            });

            requestAnimationFrame(animate);
        }

        animate();
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

            /* Prevent MathJax elements from being handled by other scripts */
            .enhanced-formula .MathJax,
            .enhanced-formula mjx-container,
            .enhanced-formula [class*="MathJax"] {
                pointer-events: none !important;
                cursor: inherit !important;
            }

            /* Ensure images inside formula wrapper are not selected by generic image handlers */
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
            @keyframes formulaZoom {
                from {
                    opacity: 0;
                    transform: scale(0.5) translateY(100px) rotateX(30deg);
                }
                to {
                    opacity: 1;
                    transform: scale(1) translateY(0) rotateX(0deg);
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

                .formula-modal {
                    padding: 1rem !important;
                }

                .formula-modal mjx-container,
                .formula-modal .MathJax {
                    font-size: 1.5rem !important;
                }
            }

            /* Dark mode support */
            @media (prefers-color-scheme: dark) {
                .enhanced-formula {
                    background: linear-gradient(135deg, #2d3748 0%, #1a202c 100%);
                    border-color: rgba(255, 255, 255, 0.1);
                    box-shadow: 0 4px 15px rgba(0, 0, 0, 0.3);
                }

                .formula-modal {
                    background: #2d3748 !important;
                    color: white !important;
                }

                .formula-modal-overlay {
                    background: rgba(0, 0, 0, 0.95) !important;
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
