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

    // ========== 中间Logo放大动画 ==========
    function addMainLogoAnimation() {
        // 查找文档主要内容区域中的logo图片
        const mainContentLogos = document.querySelectorAll('.bd-content img, .bd-article img, article img, .document img');

        mainContentLogos.forEach(img => {
            // 检查是否是logo图片（通过路径或alt文本判断）
            if (img.src && (img.src.includes('SkyBornLogo') || img.alt.toLowerCase().includes('logo'))) {
                img.style.cursor = 'pointer';

                // 移除默认的图片点击事件
                img.onclick = null;

                img.addEventListener('click', function(e) {
                    e.preventDefault();
                    e.stopPropagation();

                    createLogoModal(this);
                });
            }
        });
    }

    // 创建Logo模态弹窗
    function createLogoModal(logoElement) {
        // 防止重复创建
        if (document.querySelector('.logo-modal')) {
            return;
        }

        // 创建遮罩层
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

        // 创建模态容器
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

        // 创建关闭按钮
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

        // 创建logo图片副本
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

        // 创建描述
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

        // 组装模态内容（不包含标题）
        modal.appendChild(closeBtn);
        modal.appendChild(logoClone);
        modal.appendChild(description);
        overlay.appendChild(modal);
        document.body.appendChild(overlay);

        // 阻止页面滚动
        document.body.style.overflow = 'hidden';

        // 动画显示
        requestAnimationFrame(() => {
            overlay.style.opacity = '1';
            modal.style.transform = 'scale(1) rotate(0deg)';
            modal.style.opacity = '1';
        });

        // 关闭功能
        function closeModal() {
            overlay.style.opacity = '0';
            modal.style.transform = 'scale(0.3) rotate(-10deg)';
            modal.style.opacity = '0';

            setTimeout(() => {
                document.body.removeChild(overlay);
                document.body.style.overflow = '';
            }, 400);
        }

        // 绑定关闭事件
        closeBtn.addEventListener('click', closeModal);
        overlay.addEventListener('click', (e) => {
            if (e.target === overlay) {
                closeModal();
            }
        });

        // ESC键关闭
        const escHandler = (e) => {
            if (e.key === 'Escape') {
                closeModal();
                document.removeEventListener('keydown', escHandler);
            }
        };
        document.addEventListener('keydown', escHandler);

        // 添加粒子效果背景
        addModalParticles(modal);
    }

    // 为模态添加粒子效果
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

        // 调整canvas大小
        function resizeCanvas() {
            canvas.width = modal.offsetWidth;
            canvas.height = modal.offsetHeight;
        }
        resizeCanvas();

        // 创建粒子
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

        // 动画循环
        function animate() {
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            particles.forEach(particle => {
                particle.x += particle.vx;
                particle.y += particle.vy;

                // 边界检测
                if (particle.x < 0 || particle.x > canvas.width) particle.vx *= -1;
                if (particle.y < 0 || particle.y > canvas.height) particle.vy *= -1;

                // 绘制粒子
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
    function addLogoClickToHome() {
        // 查找所有可能的logo元素
        const logoSelectors = [
            '.bd-header .navbar-brand',
            '.bd-header .navbar-brand img',
            '.bd-header img',
            '.bd-header .logo'
        ];

        logoSelectors.forEach(selector => {
            const logos = document.querySelectorAll(selector);
            logos.forEach(logo => {
                logo.style.cursor = 'pointer';
                logo.addEventListener('click', function(e) {
                    e.preventDefault();

                    // 添加点击动画
                    this.style.transform = 'scale(0.95)';
                    this.style.transition = 'transform 0.1s ease';

                    // 重置动画
                    setTimeout(() => {
                        this.style.transform = 'scale(1)';
                    }, 100);

                    // 导航到首页（粒子效果页面）
                    setTimeout(() => {
                        window.location.href = 'index.html';
                    }, 200);
                });

                // 添加鼠标悬停效果
                logo.addEventListener('mouseenter', function() {
                    this.style.cursor = 'pointer';
                });
            });
        });
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
        addMainLogoAnimation(); // 添加中间logo动画

        console.log('🎉 Skyborn Documentation interactive effects loaded!');
    }

    // Delayed initialization to ensure all elements are loaded
    setTimeout(initializeEffects, 100);

    // 多次尝试添加logo点击事件，因为logo可能是动态加载的
    setTimeout(() => {
        addLogoClickToHome();
        addMainLogoAnimation(); // 也要重新添加中间logo事件
        console.log('Logo click handlers re-added after 1 second');
    }, 1000);

    setTimeout(() => {
        addLogoClickToHome();
        addMainLogoAnimation(); // 也要重新添加中间logo事件
        console.log('Logo click handlers re-added after 2 seconds');
    }, 2000);

    // 监听页面变化，动态添加logo点击事件
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
                addMainLogoAnimation(); // 也要重新添加中间logo事件
                console.log('Logo click handlers added after DOM mutation');
            }, 100);
        }
    });

    observer.observe(document.body, {
        childList: true,
        subtree: true
    });

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

        @keyframes logoFloat {
            0%, 100% {
                transform: translateY(0px);
            }
            50% {
                transform: translateY(-10px);
            }
        }

        /* 防止中间logo的默认点击行为 */
        .bd-content img[src*="SkyBornLogo"],
        .bd-article img[src*="SkyBornLogo"],
        article img[src*="SkyBornLogo"],
        .document img[src*="SkyBornLogo"] {
            cursor: pointer !important;
        }

        /* Enhanced button hover states */
        .btn:hover, .bd-button:hover {
            animation: pulse 2s infinite;
        }
    `;
    document.head.appendChild(style);
});
