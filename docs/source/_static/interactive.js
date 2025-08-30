/* Skyborn Documentation - Interactive Effects */

// Wait for DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {

    // ========== Mouse Follower Effect ==========
    function createMouseFollower() {
        const follower = document.createElement('div');
        follower.className = 'mouse-follower';
        document.body.appendChild(follower);

        let mouseX = 0, mouseY = 0;
        let followerX = 0, followerY = 0;

        document.addEventListener('mousemove', (e) => {
            mouseX = e.clientX;
            mouseY = e.clientY;
        });

        function updateFollower() {
            // Smooth following with easing
            followerX += (mouseX - followerX) * 0.1;
            followerY += (mouseY - followerY) * 0.1;

            follower.style.left = followerX + 'px';
            follower.style.top = followerY + 'px';

            requestAnimationFrame(updateFollower);
        }
        updateFollower();
    }

    // ========== Scroll Progress Indicator ==========
    function createScrollProgress() {
        const progressBar = document.createElement('div');
        progressBar.className = 'scroll-progress';
        document.body.appendChild(progressBar);

        window.addEventListener('scroll', () => {
            const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
            const scrollHeight = document.documentElement.scrollHeight - window.innerHeight;
            const scrollPercent = (scrollTop / scrollHeight) * 100;

            progressBar.style.width = scrollPercent + '%';
        });
    }

    // ========== Page Loading Animation ==========
    const articles = document.querySelectorAll('.bd-article');
    articles.forEach((article, index) => {
        article.style.opacity = '0';
        article.style.transform = 'translateY(40px)';

        setTimeout(() => {
            article.style.transition = 'all 0.8s cubic-bezier(0.175, 0.885, 0.32, 1.275)';
            article.style.opacity = '1';
            article.style.transform = 'translateY(0)';
        }, index * 150);
    });

    // ========== Button Click Animation ==========
    function addButtonAnimation() {
        const buttons = document.querySelectorAll('.btn, .bd-button, button');

        buttons.forEach(button => {
            button.addEventListener('click', function(e) {
                // Create ripple effect
                const ripple = document.createElement('span');
                const rect = this.getBoundingClientRect();
                const size = Math.max(rect.width, rect.height);
                const x = e.clientX - rect.left - size / 2;
                const y = e.clientY - rect.top - size / 2;

                ripple.style.width = ripple.style.height = size + 'px';
                ripple.style.left = x + 'px';
                ripple.style.top = y + 'px';
                ripple.style.position = 'absolute';
                ripple.style.borderRadius = '50%';
                ripple.style.background = 'rgba(255, 255, 255, 0.6)';
                ripple.style.transform = 'scale(0)';
                ripple.style.animation = 'ripple 0.6s linear';
                ripple.style.pointerEvents = 'none';

                this.style.position = 'relative';
                this.style.overflow = 'hidden';
                this.appendChild(ripple);

                // Remove ripple element after animation
                setTimeout(() => {
                    ripple.remove();
                }, 600);
            });
        });
    }

    // ========== Fullscreen Animation =========
    function addFullscreenAnimation() {
        const fullscreenBtn = document.querySelector('.bd-header-fullscreen');

        if (fullscreenBtn) {
            fullscreenBtn.addEventListener('click', function(e) {
                e.preventDefault();

                // Create loading state for button
                const originalText = this.innerHTML;
                this.innerHTML = '<i class="fas fa-spinner fa-spin"></i>';
                this.style.pointerEvents = 'none';

                // Create overlay for smooth transition
                const overlay = document.createElement('div');
                overlay.style.position = 'fixed';
                overlay.style.top = '0';
                overlay.style.left = '0';
                overlay.style.width = '100%';
                overlay.style.height = '100%';
                overlay.style.background = 'rgba(0, 0, 0, 0.8)';
                overlay.style.zIndex = '9998';
                overlay.style.opacity = '0';
                overlay.style.transition = 'opacity 0.3s ease';
                document.body.appendChild(overlay);

                // Animate overlay appearance
                setTimeout(() => {
                    overlay.style.opacity = '1';
                }, 10);

                // Scale animation for page content
                const mainContent = document.querySelector('.bd-container');
                if (mainContent) {
                    mainContent.style.transition = 'transform 0.4s cubic-bezier(0.25, 0.46, 0.45, 0.94)';
                    mainContent.style.transform = 'scale(0.95)';
                }

                // Simulate fullscreen toggle after animation
                setTimeout(() => {
                    // Remove overlay
                    overlay.style.opacity = '0';
                    setTimeout(() => {
                        overlay.remove();
                    }, 300);

                    // Reset button state
                    this.innerHTML = originalText;
                    this.style.pointerEvents = 'auto';

                    // Reset content scale
                    if (mainContent) {
                        mainContent.style.transform = 'scale(1)';
                    }

                    // Toggle actual fullscreen
                    if (!document.fullscreenElement) {
                        document.documentElement.requestFullscreen().catch(err => {
                            console.log('Fullscreen request failed:', err);
                        });
                    } else {
                        document.exitFullscreen();
                    }
                }, 600);
            });
        }
    }

    // ========== Enhanced Image Hover 3D Effects ==========
    function addImageEffects() {
        const images = document.querySelectorAll('.nboutput img, .cell_output img');

        images.forEach(img => {
            // Add mouse follower shadow effect
            let shadowElement = null;

            img.addEventListener('mouseenter', function() {
                this.style.transform = 'scale(1.05) rotateY(5deg)';
                this.style.transition = 'all 0.4s cubic-bezier(0.175, 0.885, 0.32, 1.275)';

                // Create dynamic shadow element
                shadowElement = document.createElement('div');
                shadowElement.style.position = 'absolute';
                shadowElement.style.width = this.offsetWidth + 'px';
                shadowElement.style.height = this.offsetHeight + 'px';
                shadowElement.style.background = 'radial-gradient(circle, rgba(37, 99, 235, 0.3) 0%, transparent 70%)';
                shadowElement.style.borderRadius = '16px';
                shadowElement.style.pointerEvents = 'none';
                shadowElement.style.zIndex = '-1';
                shadowElement.style.transition = 'all 0.1s ease';

                this.parentElement.style.position = 'relative';
                this.parentElement.appendChild(shadowElement);
            });

            img.addEventListener('mouseleave', function() {
                this.style.transform = 'scale(1) rotateY(0deg)';
                if (shadowElement) {
                    shadowElement.remove();
                    shadowElement = null;
                }
            });

            // Enhanced mouse move parallax effect with shadow following
            img.addEventListener('mousemove', function(e) {
                const rect = this.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;
                const centerX = rect.width / 2;
                const centerY = rect.height / 2;
                const rotateX = (y - centerY) / 15;
                const rotateY = (centerX - x) / 15;

                // Enhanced 3D transform
                this.style.transform = `
                    scale(1.05)
                    rotateX(${rotateX}deg)
                    rotateY(${rotateY}deg)
                    translateZ(20px)
                `;

                // Move shadow to follow mouse
                if (shadowElement) {
                    const shadowOffsetX = (x - centerX) * 0.1;
                    const shadowOffsetY = (y - centerY) * 0.1;
                    shadowElement.style.transform = `translate(${shadowOffsetX}px, ${shadowOffsetY}px)`;
                    shadowElement.style.background = `radial-gradient(circle at ${x}px ${y}px, rgba(37, 99, 235, 0.4) 0%, transparent 70%)`;
                }
            });
        });
    }

    // ========== Sidebar Navigation Animation ==========
    function addNavigationEffects() {
        const navLinks = document.querySelectorAll('.bd-sidebar-primary a');

        navLinks.forEach((link, index) => {
            // Delayed loading animation
            setTimeout(() => {
                link.style.opacity = '0';
                link.style.transform = 'translateX(-20px)';
                link.style.transition = 'all 0.3s ease-out';

                setTimeout(() => {
                    link.style.opacity = '1';
                    link.style.transform = 'translateX(0)';
                }, 50);
            }, index * 50);

            // Hover effect
            link.addEventListener('mouseenter', function() {
                this.style.transform = 'translateX(8px)';
                this.style.transition = 'all 0.3s cubic-bezier(0.175, 0.885, 0.32, 1.275)';
            });

            link.addEventListener('mouseleave', function() {
                this.style.transform = 'translateX(0)';
            });
        });
    }

    // ========== Code Block Interactive Effects ==========
    function addCodeBlockEffects() {
        const codeBlocks = document.querySelectorAll('.highlight pre');

        codeBlocks.forEach(block => {
            // Add copy button effect
            block.addEventListener('mouseenter', function() {
                this.style.transform = 'translateY(-2px)';
                this.style.boxShadow = '0 8px 25px rgba(37, 99, 235, 0.15)';
                this.style.transition = 'all 0.3s ease';
            });

            block.addEventListener('mouseleave', function() {
                this.style.transform = 'translateY(0)';
                this.style.boxShadow = '0 2px 8px rgba(0, 0, 0, 0.05)';
            });
        });
    }

    // ========== Smooth Scrolling for Anchor Links ==========
    function addSmoothScrolling() {
        const anchors = document.querySelectorAll('a[href^="#"]');

        anchors.forEach(anchor => {
            anchor.addEventListener('click', function(e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));

                if (target) {
                    target.scrollIntoView({
                        behavior: 'smooth',
                        block: 'start'
                    });

                    // Highlight target element
                    target.style.transition = 'all 0.3s ease';
                    target.style.backgroundColor = 'rgba(37, 99, 235, 0.1)';
                    target.style.padding = '10px';
                    target.style.borderRadius = '8px';

                    setTimeout(() => {
                        target.style.backgroundColor = '';
                        target.style.padding = '';
                    }, 2000);
                }
            });
        });
    }

    // ========== Scroll Progress Bar ==========
    function addScrollProgress() {
        const progressBar = document.createElement('div');
        progressBar.style.position = 'fixed';
        progressBar.style.top = '0';
        progressBar.style.left = '0';
        progressBar.style.width = '0%';
        progressBar.style.height = '3px';
        progressBar.style.background = 'linear-gradient(90deg, #2563eb, #3b82f6)';
        progressBar.style.zIndex = '9999';
        progressBar.style.transition = 'width 0.1s ease';
        document.body.appendChild(progressBar);

        window.addEventListener('scroll', function() {
            const scrollTop = window.pageYOffset;
            const docHeight = document.body.scrollHeight - window.innerHeight;
            const scrollPercent = (scrollTop / docHeight) * 100;
            progressBar.style.width = scrollPercent + '%';
        });
    }

    // ========== Theme Transition Animations ==========
    function addThemeTransitions() {
        const style = document.createElement('style');
        style.textContent = `
            @keyframes ripple {
                to {
                    transform: scale(4);
                    opacity: 0;
                }
            }

            @keyframes float {
                0%, 100% { transform: translateY(0px); }
                50% { transform: translateY(-10px); }
            }

            .float-animation {
                animation: float 3s ease-in-out infinite;
            }
        `;
        document.head.appendChild(style);
    }

    // ========== Image Zoom Utility ==========
    function addImageZoomFunction() {
    // Find all images (including Jupyter notebook outputs and documentation images)
    // but exclude MathJax-related images and elements
        const allImages = document.querySelectorAll('img:not(.MathJax img):not([class*="MathJax"]):not(mjx-container img)');

        allImages.forEach(img => {
            // Skip navbar logos, favicons, MathJax-related small images, etc.
            if (img.src &&
                !img.src.includes('favicon') &&
                !img.classList.contains('navbar-logo') &&
                !img.closest('.MathJax') && // exclude images inside MathJax containers
                !img.closest('mjx-container') && // exclude images inside mjx-container
                !img.closest('.enhanced-formula') && // exclude images inside enhanced-formula wrappers
                !img.parentElement.classList.contains('formula-wrapper') && // exclude formula wrapper parents
                img.offsetWidth > 100 && // only handle larger images
                img.offsetHeight > 100) {

                img.style.cursor = 'pointer';
                img.style.transition = 'transform 0.3s ease';

                // add hover effect
                img.addEventListener('mouseenter', function() {
                    this.style.transform = 'scale(1.05)';
                });

                img.addEventListener('mouseleave', function() {
                    this.style.transform = 'scale(1)';
                });

                // remove default image click handler
                img.onclick = null;

                img.addEventListener('click', function(e) {
                    e.preventDefault();
                    e.stopPropagation();

                    // determine if it's a logo and use a different modal
                    if (this.src.includes('SkyBornLogo') || this.alt.toLowerCase().includes('logo')) {
                        createLogoModal(this);
                    } else {
                        createImageModal(this);
                    }
                });
            }
        });
    }

    // ========== Center Logo Zoom Animation ==========
    function addMainLogoAnimation() {
    // This function is now handled by addImageZoomFunction
    // Kept for backward compatibility
    }

    // Create logo modal dialog
    function createLogoModal(logoElement) {
    // prevent duplicate creation
        if (document.querySelector('.logo-modal')) {
            return;
        }

    // create overlay
        const overlay = document.createElement('div');
        overlay.className = 'logo-modal-overlay';
        overlay.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0, 0, 0, 0.8);
            z-index: 10000;
            opacity: 0;
            transition: opacity 0.3s ease;
            display: flex;
            align-items: center;
            justify-content: center;
            backdrop-filter: blur(10px);
        `;

    // create modal container
        const modal = document.createElement('div');
        modal.className = 'logo-modal';
        modal.style.cssText = `
            background: white;
            border-radius: 20px;
            padding: 2rem;
            width: 50vw;
            height: 50vh;
            min-width: 400px;
            min-height: 300px;
            max-width: 800px;
            max-height: 600px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            transform: scale(0.3) rotate(10deg);
            transition: all 0.4s cubic-bezier(0.175, 0.885, 0.32, 1.275);
            opacity: 0;
            position: relative;
            overflow: hidden;
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
        `;

    // create close button
        const closeBtn = document.createElement('button');
        closeBtn.innerHTML = '×';
        closeBtn.style.cssText = `
            position: absolute;
            top: 10px;
            right: 15px;
            background: none;
            border: none;
            font-size: 2rem;
            cursor: pointer;
            color: #999;
            transition: color 0.3s ease;
            z-index: 10001;
        `;

        closeBtn.addEventListener('mouseenter', () => {
            closeBtn.style.color = '#333';
            closeBtn.style.transform = 'scale(1.2)';
        });

        closeBtn.addEventListener('mouseleave', () => {
            closeBtn.style.color = '#999';
            closeBtn.style.transform = 'scale(1)';
        });

    // create a clone of the logo image
        const logoClone = logoElement.cloneNode(true);
        logoClone.style.cssText = `
            max-width: 70%;
            max-height: 70%;
            border: none !important;
            border-radius: 0 !important;
            padding: 0 !important;
            background: none !important;
            box-shadow: none !important;
            margin: 0 auto;
            display: block;
            animation: logoFloat 3s ease-in-out infinite;
        `;

    // create description
        const description = document.createElement('p');
        description.textContent = 'Climate Data Analysis & Visualization Toolkit';
        description.style.cssText = `
            text-align: center;
            color: #666;
            margin-top: 2rem;
            font-style: italic;
            font-size: 1.2rem;
            font-weight: 500;
        `;

    // assemble modal content (excluding title)
        modal.appendChild(closeBtn);
        modal.appendChild(logoClone);
        modal.appendChild(description);
        overlay.appendChild(modal);
        document.body.appendChild(overlay);

    // prevent page scrolling
        document.body.style.overflow = 'hidden';

    // animate display
        requestAnimationFrame(() => {
            overlay.style.opacity = '1';
            modal.style.transform = 'scale(1) rotate(0deg)';
            modal.style.opacity = '1';
        });

    // close functionality
        function closeModal() {
            overlay.style.opacity = '0';
            modal.style.transform = 'scale(0.3) rotate(-10deg)';
            modal.style.opacity = '0';

            setTimeout(() => {
                document.body.removeChild(overlay);
                document.body.style.overflow = '';
            }, 400);
        }

    // bind close events
        closeBtn.addEventListener('click', closeModal);
        overlay.addEventListener('click', (e) => {
            if (e.target === overlay) {
                closeModal();
            }
        });

    // close on ESC key
        const escHandler = (e) => {
            if (e.key === 'Escape') {
                closeModal();
                document.removeEventListener('keydown', escHandler);
            }
        };
        document.addEventListener('keydown', escHandler);

    // add particle effect background
        addModalParticles(modal);
    }

    // Add modal particle effects
    function addModalParticles(modal) {
        const canvas = document.createElement('canvas');
        canvas.style.cssText = `
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            pointer-events: none;
            z-index: -1;
        `;

        modal.appendChild(canvas);

        const ctx = canvas.getContext('2d');
        const particles = [];

    // resize canvas
        function resizeCanvas() {
            canvas.width = modal.offsetWidth;
            canvas.height = modal.offsetHeight;
        }
        resizeCanvas();

    // create particles
        for (let i = 0; i < 30; i++) {
            particles.push({
                x: Math.random() * canvas.width,
                y: Math.random() * canvas.height,
                vx: (Math.random() - 0.5) * 2,
                vy: (Math.random() - 0.5) * 2,
                size: Math.random() * 3 + 1,
                opacity: Math.random() * 0.5 + 0.3,
                color: `hsl(${210 + Math.random() * 30}, 70%, 60%)`
            });
        }

    // animation loop
        function animate() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            particles.forEach(particle => {
                particle.x += particle.vx;
                particle.y += particle.vy;

                // boundary detection
                if (particle.x < 0 || particle.x > canvas.width) particle.vx *= -1;
                if (particle.y < 0 || particle.y > canvas.height) particle.vy *= -1;

                // draw particles
                ctx.beginPath();
                ctx.arc(particle.x, particle.y, particle.size, 0, Math.PI * 2);
                ctx.fillStyle = particle.color;
                ctx.globalAlpha = particle.opacity;
                ctx.fill();
            });

            if (modal.parentNode) {
                requestAnimationFrame(animate);
            }
        }
        animate();
    }

    // ========== Generic Image Modal Dialog ==========
    function createImageModal(imageElement) {
    // prevent duplicate creation
        if (document.querySelector('.image-modal')) {
            return;
        }

    // create overlay
        const overlay = document.createElement('div');
        overlay.className = 'image-modal-overlay';
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

    // create modal container
        const modal = document.createElement('div');
        modal.className = 'image-modal';
        modal.style.cssText = `
            background: white;
            border-radius: 20px;
            padding: 1rem;
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

    // create close button
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

    // create image clone
        const imageClone = imageElement.cloneNode(true);
        imageClone.style.cssText = `
            max-width: 100%;
            max-height: 100%;
            width: auto;
            height: auto;
            border: none !important;
            border-radius: 10px !important;
            padding: 0 !important;
            background: none !important;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.3) !important;
            margin: 0 auto;
            display: block;
            object-fit: contain;
            animation: imageZoom 0.5s ease-out;
        `;

    // create image info
        const imageInfo = document.createElement('div');
        const fileName = imageElement.src.split('/').pop() || 'Image';
        const imageTitle = imageElement.alt || imageElement.title || fileName;

        // imageInfo.innerHTML = `
        //     <h3 style="margin: 1rem 0 0.5rem; color: #333; font-size: 1.2rem; text-align: center;">${imageTitle}</h3>
        //     <p style="margin: 0; color: #666; font-size: 0.9rem; text-align: center;">Click and drag to pan • Scroll to zoom</p>
        // `;

    // add image interaction features
        let isDragging = false;
        let startX, startY, translateX = 0, translateY = 0;
        let scale = 1;
        let lastTranslateX = 0, lastTranslateY = 0;

    // mouse down: start dragging
        imageClone.addEventListener('mousedown', (e) => {
            e.preventDefault();
            isDragging = true;
            startX = e.clientX - translateX;
            startY = e.clientY - translateY;
            imageClone.style.cursor = 'grabbing';
            imageClone.style.userSelect = 'none';
            imageClone.style.transition = 'none'; // disable transitions during drag
        });

    // mouse move: dragging image
        document.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            e.preventDefault();
            translateX = e.clientX - startX;
            translateY = e.clientY - startY;
            imageClone.style.transform = `translate(${translateX}px, ${translateY}px) scale(${scale})`;
        });

    // mouse up: end dragging
        document.addEventListener('mouseup', (e) => {
            if (isDragging) {
                isDragging = false;
                imageClone.style.cursor = 'grab';
                imageClone.style.userSelect = 'auto';
                imageClone.style.transition = 'transform 0.2s ease'; // restore transitions
                // save current position
                lastTranslateX = translateX;
                lastTranslateY = translateY;
            }
        });

    // wheel zoom (optimize zoom center)
        imageClone.addEventListener('wheel', (e) => {
            e.preventDefault();
            const rect = imageClone.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            const delta = e.deltaY > 0 ? 0.9 : 1.1;
            const newScale = Math.max(0.5, Math.min(5, scale * delta));

            if (newScale !== scale) {
                // calculate zoom center adjustment
                const scaleChange = newScale / scale;
                translateX = translateX * scaleChange + mouseX * (1 - scaleChange);
                translateY = translateY * scaleChange + mouseY * (1 - scaleChange);

                scale = newScale;
                imageClone.style.transform = `translate(${translateX}px, ${translateY}px) scale(${scale})`;

                // update saved position
                lastTranslateX = translateX;
                lastTranslateY = translateY;
            }
        });

    // double-click: reset position and zoom
        imageClone.addEventListener('dblclick', (e) => {
            e.preventDefault();
            translateX = 0;
            translateY = 0;
            scale = 1;
            imageClone.style.transition = 'transform 0.3s ease';
            imageClone.style.transform = `translate(0px, 0px) scale(1)`;
            setTimeout(() => {
                imageClone.style.transition = 'transform 0.2s ease';
            }, 300);
            lastTranslateX = 0;
            lastTranslateY = 0;
        });

    // prevent default drag behavior for images
        imageClone.addEventListener('dragstart', (e) => {
            e.preventDefault();
        });

        imageClone.style.cursor = 'grab';

    // assemble modal content
        modal.appendChild(closeBtn);
        modal.appendChild(imageClone);
        modal.appendChild(imageInfo);
        overlay.appendChild(modal);
        document.body.appendChild(overlay);

    // prevent page scrolling
        document.body.style.overflow = 'hidden';

    // animate display
        requestAnimationFrame(() => {
            overlay.style.opacity = '1';
            modal.style.transform = 'scale(1) rotate(0deg)';
            modal.style.opacity = '1';
        });

    // close functionality
        function closeModal() {
            overlay.style.opacity = '0';
            modal.style.transform = 'scale(0.3) rotate(-5deg)';
            modal.style.opacity = '0';

            setTimeout(() => {
                document.body.removeChild(overlay);
                document.body.style.overflow = '';
            }, 400);
        }

    // bind close events
        closeBtn.addEventListener('click', closeModal);
        overlay.addEventListener('click', (e) => {
            if (e.target === overlay) {
                closeModal();
            }
        });

    // close on ESC key
        const escHandler = (e) => {
            if (e.key === 'Escape') {
                closeModal();
                document.removeEventListener('keydown', escHandler);
            }
        };
        document.addEventListener('keydown', escHandler);

    // add background animation effects
        addImageModalEffects(modal);
    }

    // Add background effects for image modal
    function addImageModalEffects(modal) {
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

    // create decorative dots
        for (let i = 0; i < 20; i++) {
            dots.push({
                x: Math.random() * canvas.width,
                y: Math.random() * canvas.height,
                vx: (Math.random() - 0.5) * 1,
                vy: (Math.random() - 0.5) * 1,
                size: Math.random() * 2 + 1,
                opacity: Math.random() * 0.3 + 0.1,
                color: `hsl(${200 + Math.random() * 60}, 70%, 70%)`
            });
        }

        function animate() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            dots.forEach(dot => {
                dot.x += dot.vx;
                dot.y += dot.vy;

                if (dot.x < 0 || dot.x > canvas.width) dot.vx *= -1;
                if (dot.y < 0 || dot.y > canvas.height) dot.vy *= -1;

                ctx.beginPath();
                ctx.arc(dot.x, dot.y, dot.size, 0, Math.PI * 2);
                ctx.fillStyle = dot.color;
                ctx.globalAlpha = dot.opacity;
                ctx.fill();
            });

            if (modal.parentNode) {
                requestAnimationFrame(animate);
            }
        }
        animate();
    }
    function addLogoClickToHome() {
    // find all possible logo elements
        const logoSelectors = [
            '.bd-header .navbar-brand',
            '.bd-header .navbar-brand img',
            '.bd-header img',
            '.bd-header .logo',
            '.bd-sidebar-primary .navbar-brand',
            '.bd-sidebar-primary .navbar-brand img',
            '.navbar-brand',
            '.navbar-brand img',
            'a.navbar-brand.logo',
            'a.navbar-brand.logo img'
        ];

        logoSelectors.forEach(selector => {
            const logos = document.querySelectorAll(selector);
            logos.forEach(logo => {
                // set styles
                logo.style.cursor = 'pointer';

                // remove previous event listeners (prevent duplicate bindings)
                logo.removeEventListener('click', handleLogoClick);

                // add new event listeners
                logo.addEventListener('click', handleLogoClick);

                // add hover effect
                logo.addEventListener('mouseenter', function() {
                    this.style.cursor = 'pointer';
                });

                // force pointer style
                logo.classList.add('logo-clickable');
            });
        });

    // Logo click handler function
        function handleLogoClick(e) {
            e.preventDefault();
            e.stopPropagation();

            console.log('Logo clicked! Navigating to index.html');

            // add click animation effect
            this.style.transform = 'scale(0.9)';
            this.style.transition = 'transform 0.15s ease';

            // add ripple effect
            const ripple = document.createElement('div');
            const rect = this.getBoundingClientRect();
            const size = Math.min(rect.width, rect.height) * 0.8; // scaled down to 0.8x

            ripple.style.cssText = `
                position: fixed;
                left: ${rect.left + rect.width/2 - size/2}px;
                top: ${rect.top + rect.height/2 - size/2}px;
                width: ${size}px;
                height: ${size}px;
                border-radius: 50%;
                background: radial-gradient(circle, rgba(59, 130, 246, 0.4) 0%, rgba(147, 197, 253, 0.2) 50%, rgba(59, 130, 246, 0) 100%);
                transform: scale(0);
                animation: ripple 0.5s ease-out;
                pointer-events: none;
                z-index: 10000;
            `;

            document.body.appendChild(ripple);

            // reset after animation completes
            setTimeout(() => {
                this.style.transform = 'scale(1)';
                if (ripple.parentNode) {
                    ripple.remove();
                }
            }, 150);

            // navigate to homepage (particle effect page)
            setTimeout(() => {
                window.location.href = 'index.html';
            }, 200);
        }

        console.log('🎯 Logo click handlers added successfully');
    }

    // ========== Initialize All Effects ==========
    function initializeEffects() {
        addButtonAnimation();
        addFullscreenAnimation();
        addImageEffects();
        addNavigationEffects();
        addCodeBlockEffects();
        addSmoothScrolling();
        addScrollProgress();
        addThemeTransitions();
        addLogoClickToHome();
    addImageZoomFunction(); // add generic image zoom function (including logos)

        console.log('🎉 Skyborn Documentation interactive effects loaded!');
    }

    // Call new functions
    createMouseFollower();
    createScrollProgress();

    // Delayed initialization to ensure all elements are loaded
    setTimeout(initializeEffects, 100);

    // Retry adding logo click handlers multiple times because logos may load dynamically
    setTimeout(() => {
        addLogoClickToHome();
    addImageZoomFunction(); // also re-add image zoom handlers
        console.log('Logo click handlers re-added after 1 second');
    }, 1000);

    setTimeout(() => {
        addLogoClickToHome();
    addImageZoomFunction(); // also re-add image zoom handlers
        console.log('Logo click handlers re-added after 2 seconds');
    }, 2000);

    // Observe DOM mutations to dynamically add logo click handlers
    const observer = new MutationObserver(function(mutations) {
        let shouldAddLogoHandler = false;
        mutations.forEach(function(mutation) {
            if (mutation.type === 'childList' && mutation.addedNodes.length > 0) {
                mutation.addedNodes.forEach(function(node) {
                    if (node.nodeType === Node.ELEMENT_NODE) {
                        if (node.classList && (node.classList.contains('bd-header') ||
                            node.classList.contains('navbar-brand') ||
                            node.tagName === 'IMG')) {
                            shouldAddLogoHandler = true;
                        }
                    }
                });
            }
        });

        if (shouldAddLogoHandler) {
            setTimeout(() => {
                addLogoClickToHome();
                addImageZoomFunction(); // also re-add image zoom handlers
                console.log('Logo click handlers added after DOM mutation');
            }, 100);
        }
    });

    observer.observe(document.body, {
        childList: true,
        subtree: true
    });


    // ========== Additional CSS Animations ==========
    const style = document.createElement('style');
    style.textContent = `
        @keyframes ripple {
            to {
                transform: scale(4);
                opacity: 0;
            }
        }

        @keyframes fadeInUp {
            from {
                opacity: 0;
                transform: translateY(40px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        @keyframes pulse {
            0%, 100% {
                transform: scale(1);
            }
            50% {
                transform: scale(1.05);
            }
        }

        @keyframes float {
            0%, 100% {
                transform: translateY(0px);
            }
            50% {
                transform: translateY(-10px);
            }
        }

        /* Perspective for enhanced 3D effects */
        .nboutput, .cell_output {
            perspective: 1200px;
            transform-style: preserve-3d;
        }

        .nboutput img, .cell_output img {
            transform-style: preserve-3d;
        }

        @keyframes logoFloat {
            0%, 100% {
                transform: translateY(0px);
            }
            50% {
                transform: translateY(-10px);
            }
        }

        @keyframes imageZoom {
            0% {
                transform: scale(0.8) rotate(5deg);
                opacity: 0;
            }
            100% {
                transform: scale(1) rotate(0deg);
                opacity: 1;
            }
        }

    /* General image hover effects */
        img:not(.navbar-brand img):not([src*="favicon"]) {
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }

        img:not(.navbar-brand img):not([src*="favicon"]):hover {
            box-shadow: 0 10px 25px rgba(0, 0, 0, 0.2);
        }

    /* Prevent default click behavior for center logo */
        .bd-content img[src*="SkyBornLogo"],
        .bd-article img[src*="SkyBornLogo"],
        article img[src*="SkyBornLogo"],
        .document img[src*="SkyBornLogo"] {
            cursor: pointer !important;
        }

    /* All large images are clickable for zoom */
        img[width]:not(.navbar-brand img):not([src*="favicon"]) {
            cursor: pointer !important;
        }

    /* Styles for image modal window */
        .image-modal-overlay {
            backdrop-filter: blur(15px) !important;
        }

        .image-modal {
            background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%) !important;
            border: 1px solid rgba(255, 255, 255, 0.2);
        }

        /* Enhanced button hover states */
        .btn:hover, .bd-button:hover {
            animation: pulse 2s infinite;
        }
    `;
    document.head.appendChild(style);
});
