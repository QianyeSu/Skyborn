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

    // ========== Initialize All Effects ==========
    function initializeEffects() {
        addButtonAnimation();
        addImageEffects();
        addNavigationEffects();
        addCodeBlockEffects();
        addSmoothScrolling();
        addScrollProgress();
        addThemeTransitions();

        console.log('🎉 Skyborn Documentation interactive effects loaded!');
    }

    // Delayed initialization to ensure all elements are loaded
    setTimeout(initializeEffects, 100);

    // Call new functions
    createMouseFollower();
    createScrollProgress();

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

        /* Enhanced button hover states */
        .btn:hover, .bd-button:hover {
            animation: pulse 2s infinite;
        }
    `;
    document.head.appendChild(style);
});
