/* Skyborn Documentation - Dynamic Table Column Width Adjustment */

document.addEventListener('DOMContentLoaded', function() {

    // ========== Table Auto-Adjustment System ==========
    function initializeTableAdjustment() {

        // 查找所有表格
        const tables = document.querySelectorAll('table');

        tables.forEach(table => {
            makeTableResponsive(table);
        });

        // 监听窗口大小变化 - 但不自动调整已手动调整的表格
        let resizeTimer;
        window.addEventListener('resize', function() {
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(() => {
                tables.forEach(table => {
                    // 只调整未被用户手动调整的表格
                    if (!table.dataset.userAdjusted) {
                        adjustTableColumns(table);
                    }
                });
            }, 250);
        });

        console.log('🔧 Table responsive system initialized');
    }

    // ========== Make Table Responsive ==========
    function makeTableResponsive(table) {

        // 跳过已经处理过的表格
        if (table.classList.contains('responsive-table-processed')) {
            return;
        }

        // 添加响应式容器
        const wrapper = document.createElement('div');
        wrapper.className = 'responsive-table-wrapper';
        wrapper.style.cssText = `
            position: relative;
            width: 100%;
            overflow: hidden;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(37, 99, 235, 0.15);
            margin: 1rem 0;
            background: white;
            transition: all 0.3s ease;
        `;

        // 创建表格容器
        const container = document.createElement('div');
        container.className = 'table-scroll-container';
        container.style.cssText = `
            width: 100%;
            overflow-x: auto;
            overflow-y: hidden;
            position: relative;
            scroll-behavior: smooth;
            /* 自定义滚动条样式 */
            scrollbar-width: thin;
            scrollbar-color: #3b82f6 #f1f5f9;
        `;

        // WebKit 浏览器滚动条样式
        const scrollbarStyle = document.createElement('style');
        scrollbarStyle.textContent = `
            .table-scroll-container::-webkit-scrollbar {
                height: 8px;
            }
            .table-scroll-container::-webkit-scrollbar-track {
                background: #f1f5f9;
                border-radius: 4px;
            }
            .table-scroll-container::-webkit-scrollbar-thumb {
                background: #3b82f6;
                border-radius: 4px;
                transition: background 0.3s ease;
            }
            .table-scroll-container::-webkit-scrollbar-thumb:hover {
                background: #2563eb;
            }
        `;
        document.head.appendChild(scrollbarStyle);

        // 包装表格
        table.parentNode.insertBefore(wrapper, table);
        wrapper.appendChild(container);
        container.appendChild(table);

        // 添加表格控制栏
        addTableControls(wrapper, table, container);

        // 初始调整表格 - 但尊重 RST 定义的列宽
        applyInitialLayout(table);

        // 添加交互效果
        addTableInteraction(table, container);

        // 标记为已处理
        table.classList.add('responsive-table-processed');

        console.log('📊 Table made responsive:', table);
    }

    // ========== Add Table Controls ==========
    function addTableControls(wrapper, table, container) {
        const controlBar = document.createElement('div');
        controlBar.className = 'table-controls';
        controlBar.style.cssText = `
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 0.5rem 1rem;
            background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%);
            border-bottom: 1px solid #e2e8f0;
            font-size: 0.85rem;
            color: #1e40af;
        `;

        // 左侧：表格标题（从前面的标题中获取）
        const info = document.createElement('div');
        const tableTitleData = getTableTitle(table);

        // 创建可点击的标题
        const titleElement = document.createElement('div');
        titleElement.innerHTML = tableTitleData.title;
        titleElement.style.cssText = `
            font-weight: 700;
            cursor: pointer;
            color: #1e40af;
            font-size: ${getHeadingSize(tableTitleData.headingLevel)};
            margin: 0;
            padding: 0.25rem 0;
            transition: color 0.3s ease;
            border-bottom: 2px solid transparent;
        `;

        // 点击跳转到原标题位置
        titleElement.addEventListener('click', () => {
            if (tableTitleData.headingElement) {
                tableTitleData.headingElement.scrollIntoView({ behavior: 'smooth', block: 'start' });
                // 高亮原标题
                tableTitleData.headingElement.style.transition = 'background-color 0.5s ease';
                tableTitleData.headingElement.style.backgroundColor = 'rgba(59, 130, 246, 0.1)';
                setTimeout(() => {
                    tableTitleData.headingElement.style.backgroundColor = '';
                }, 2000);
            }
        });

        titleElement.addEventListener('mouseenter', () => {
            titleElement.style.color = '#3b82f6';
            titleElement.style.borderBottomColor = '#3b82f6';
        });

        titleElement.addEventListener('mouseleave', () => {
            titleElement.style.color = '#1e40af';
            titleElement.style.borderBottomColor = 'transparent';
        });

        info.appendChild(titleElement);

        // 隐藏原来的标题
        if (tableTitleData.headingElement) {
            tableTitleData.headingElement.style.display = 'none';
        }

        // 右侧：控制按钮
        const controls = document.createElement('div');
        controls.style.cssText = `
            display: flex;
            gap: 0.5rem;
            align-items: center;
        `;

        // 组合控制按钮（减少按钮数量）
        const controlBtn = createControlButton('⚙️ Optimize', () => {
            // 智能优化：自动调整列宽并应用最佳设置
            autoFitColumns(table);
            table.dataset.userAdjusted = 'true'; // 标记为用户调整
            showToast('✨ Table optimized');
        });

        controls.appendChild(controlBtn);

        controlBar.appendChild(info);
        controlBar.appendChild(controls);

        wrapper.insertBefore(controlBar, container);
    }

    // ========== Get Heading Size ==========
    function getHeadingSize(headingLevel) {
        const sizes = {
            'H1': '2rem',
            'H2': '1.75rem',
            'H3': '1.5rem',
            'H4': '1.25rem',
            'H5': '1.125rem',
            'H6': '1rem'
        };
        return sizes[headingLevel] || '1.25rem';
    }

    // ========== Create Control Button ==========
    function createControlButton(text, onClick) {
        const btn = document.createElement('button');
        btn.textContent = text;
        btn.style.cssText = `
            background: rgba(255, 255, 255, 0.8);
            border: 1px solid #3b82f6;
            border-radius: 6px;
            padding: 0.25rem 0.5rem;
            font-size: 0.75rem;
            color: #1e40af;
            cursor: pointer;
            transition: all 0.2s ease;
            white-space: nowrap;
        `;

        btn.addEventListener('mouseenter', () => {
            btn.style.background = '#3b82f6';
            btn.style.color = 'white';
            btn.style.transform = 'translateY(-1px)';
            btn.style.boxShadow = '0 2px 8px rgba(59, 130, 246, 0.3)';
        });

        btn.addEventListener('mouseleave', () => {
            btn.style.background = 'rgba(255, 255, 255, 0.8)';
            btn.style.color = '#1e40af';
            btn.style.transform = 'translateY(0)';
            btn.style.boxShadow = 'none';
        });

        btn.addEventListener('click', onClick);

        return btn;
    }

    // ========== Get Table Title ==========
    function getTableTitle(table) {
        // 查找表格前面的标题
        let currentElement = table.parentElement;
        let targetHeading = null;

        // 向上查找容器中的标题
        while (currentElement && currentElement !== document.body) {
            // 查找前面的标题元素
            const prevElements = [];
            let prev = currentElement.previousElementSibling;

            // 收集前面的几个元素
            while (prev && prevElements.length < 5) {
                prevElements.unshift(prev);
                prev = prev.previousElementSibling;
            }

            // 在前面的元素中查找标题
            for (const element of prevElements) {
                // 检查是否是标题元素
                if (element.tagName && /^H[1-6]$/.test(element.tagName)) {
                    targetHeading = element;
                    return {
                        title: element.textContent.trim(),
                        headingElement: element,
                        headingLevel: element.tagName
                    };
                }

                // 检查是否包含标题文本
                if (element.classList.contains('section-header') ||
                    element.classList.contains('table-title')) {
                    return {
                        title: element.textContent.trim(),
                        headingElement: element,
                        headingLevel: 'H3'
                    };
                }

                // 查找子元素中的标题
                const heading = element.querySelector('h1, h2, h3, h4, h5, h6');
                if (heading) {
                    return {
                        title: heading.textContent.trim(),
                        headingElement: heading,
                        headingLevel: heading.tagName
                    };
                }
            }

            currentElement = currentElement.parentElement;
        }

        // 如果找不到标题，根据表格内容推断
        const firstHeaderCell = table.querySelector('th');
        if (firstHeaderCell) {
            const headerText = firstHeaderCell.textContent.trim();

            // 根据第一列标题推断表格类型
            if (headerText.includes('Function')) {
                return { title: 'Functions Reference', headingElement: null, headingLevel: 'H3' };
            } else if (headerText.includes('Class')) {
                return { title: 'Classes Reference', headingElement: null, headingLevel: 'H3' };
            } else if (headerText.includes('Method')) {
                return { title: 'Methods Reference', headingElement: null, headingLevel: 'H3' };
            } else if (headerText.includes('Parameter')) {
                return { title: 'Parameters', headingElement: null, headingLevel: 'H3' };
            } else if (headerText.includes('Attribute')) {
                return { title: 'Attributes', headingElement: null, headingLevel: 'H3' };
            } else if (headerText.includes('Page') || headerText.includes('Section')) {
                return { title: 'Documentation Structure', headingElement: null, headingLevel: 'H3' };
            }
        }

        // 检查表格是否在特定的模块section中
        const section = table.closest('[id*="calc"], [id*="spharm"], [id*="gridfill"], [id*="gradients"], [id*="causality"], [id*="windspharm"]');
        if (section) {
            const sectionId = section.id;
            if (sectionId.includes('calc')) return { title: 'Calculation Functions', headingElement: null, headingLevel: 'H3' };
            if (sectionId.includes('spharm')) return { title: 'Spherical Harmonics Functions', headingElement: null, headingLevel: 'H3' };
            if (sectionId.includes('gridfill')) return { title: 'Grid Filling Functions', headingElement: null, headingLevel: 'H3' };
            if (sectionId.includes('gradients')) return { title: 'Gradient Functions', headingElement: null, headingLevel: 'H3' };
            if (sectionId.includes('causality')) return { title: 'Causality Functions', headingElement: null, headingLevel: 'H3' };
            if (sectionId.includes('windspharm')) return { title: 'Wind Spherical Harmonics', headingElement: null, headingLevel: 'H3' };
        }

        // 默认标题
        return { title: 'Data Table', headingElement: null, headingLevel: 'H3' };
    }

    // ========== Apply Initial Layout ==========
    function applyInitialLayout(table) {
        const rows = table.querySelectorAll('tr');
        if (rows.length === 0) return;

        const firstRow = rows[0];
        const cells = firstRow.querySelectorAll('th, td');
        const columnCount = cells.length;

        if (columnCount === 0) return;

        // 为 functions_classes 表格应用预设的 30:70 比例
        if (table.closest('.skyborn-function-table') ||
            table.querySelector('th')?.textContent.includes('Function')) {

            cells.forEach((cell, index) => {
                if (index === 0) {
                    cell.style.width = '30%';
                    cell.style.maxWidth = '30%';
                } else if (index === 1) {
                    cell.style.width = '70%';
                    cell.style.maxWidth = '70%';
                }
            });

            // 应用到所有行
            rows.forEach(row => {
                const rowCells = row.querySelectorAll('th, td');
                rowCells.forEach((cell, index) => {
                    if (index === 0) {
                        cell.style.width = '30%';
                        cell.style.maxWidth = '30%';
                    } else if (index === 1) {
                        cell.style.width = '70%';
                        cell.style.maxWidth = '70%';
                    }
                });
            });

            console.log('📐 Applied 30:70 layout to function table');
        } else {
            // 对其他表格应用轻微的自动调整
            adjustTableColumns(table);
        }
    }

    // ========== Adjust Table Columns ==========
    function adjustTableColumns(table) {
        const rows = table.querySelectorAll('tr');
        if (rows.length === 0) return;

        // 获取第一行来确定列数
        const firstRow = rows[0];
        const cells = firstRow.querySelectorAll('th, td');
        const columnCount = cells.length;

        if (columnCount === 0) return;

        // 计算每列的最大内容宽度
        const columnWidths = new Array(columnCount).fill(0);
        const columnContents = new Array(columnCount).fill().map(() => []);

        rows.forEach(row => {
            const cells = row.querySelectorAll('th, td');
            cells.forEach((cell, index) => {
                if (index < columnCount) {
                    const textLength = cell.textContent.trim().length;
                    const hasLongFunction = cell.textContent.includes('skyborn.') && textLength > 25;

                    // 为长函数名分配更多空间
                    const adjustedLength = hasLongFunction ? textLength * 1.2 : textLength;
                    columnWidths[index] = Math.max(columnWidths[index], adjustedLength);
                    columnContents[index].push(cell.textContent.trim());
                }
            });
        });

        // 计算容器可用宽度
        const container = table.closest('.table-scroll-container');
        const availableWidth = container ? container.clientWidth - 40 : window.innerWidth - 100;

        // 智能分配列宽
        const totalContentWidth = columnWidths.reduce((sum, width) => sum + width, 0);
        const avgCharWidth = 8; // 估算的字符宽度
        const estimatedTableWidth = totalContentWidth * avgCharWidth;

        // 如果内容宽度超过可用空间，使用比例分配
        if (estimatedTableWidth > availableWidth) {
            const scale = availableWidth / estimatedTableWidth;
            columnWidths.forEach((width, index) => {
                columnWidths[index] = Math.max(width * scale, 100); // 最小宽度100px
            });
        }

        // 应用列宽
        applyColumnWidths(table, columnWidths);

        // 添加表格样式增强
        enhanceTableAppearance(table);
    }

    // ========== Apply Column Widths ==========
    function applyColumnWidths(table, widths) {
        // 设置表格样式
        table.style.cssText = `
            width: 100%;
            table-layout: fixed;
            border-collapse: collapse;
            margin: 0;
            background: white;
            transition: all 0.3s ease;
        `;

        // 创建或更新 colgroup
        let colgroup = table.querySelector('colgroup');
        if (!colgroup) {
            colgroup = document.createElement('colgroup');
            table.insertBefore(colgroup, table.firstChild);
        }

        colgroup.innerHTML = '';
        widths.forEach(width => {
            const col = document.createElement('col');
            col.style.width = width + 'px';
            colgroup.appendChild(col);
        });

        // 应用单元格样式
        const cells = table.querySelectorAll('th, td');
        cells.forEach(cell => {
            cell.style.cssText = `
                padding: 12px 16px;
                border: 1px solid #e2e8f0;
                text-align: left;
                vertical-align: top;
                word-wrap: break-word;
                overflow-wrap: break-word;
                hyphens: auto;
                line-height: 1.5;
                transition: all 0.2s ease;
            `;
        });

        // 头部样式
        const headers = table.querySelectorAll('th');
        headers.forEach(header => {
            header.style.cssText += `
                background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%);
                font-weight: 600;
                color: #1e40af;
                position: sticky;
                top: 0;
                z-index: 10;
            `;
        });
    }

    // ========== Auto Fit Columns ==========
    function autoFitColumns(table) {
        // 临时移除 table-layout: fixed 来测量内容
        table.style.tableLayout = 'auto';
        table.style.width = 'auto';

        // 强制浏览器重新计算
        table.offsetHeight;

        setTimeout(() => {
            // 获取自然宽度
            const rows = table.querySelectorAll('tr');
            const firstRow = rows[0];
            const cells = firstRow.querySelectorAll('th, td');
            const naturalWidths = [];

            cells.forEach(cell => {
                naturalWidths.push(cell.offsetWidth);
            });

            // 重新应用固定布局
            applyColumnWidths(table, naturalWidths);

            // 添加动画效果
            table.style.transform = 'scale(1.01)';
            setTimeout(() => {
                table.style.transform = 'scale(1)';
            }, 200);

        }, 50);
    }

    // ========== Reset Table Columns ==========
    function resetTableColumns(table) {
        // 移除自定义样式
        table.style.cssText = '';

        const colgroup = table.querySelector('colgroup');
        if (colgroup) {
            colgroup.remove();
        }

        // 重新初始化
        setTimeout(() => {
            adjustTableColumns(table);
        }, 50);
    }

    // ========== Toggle Compact Mode ==========
    function toggleCompactMode(table) {
        const isCompact = table.classList.contains('compact-mode');

        if (isCompact) {
            table.classList.remove('compact-mode');
            const cells = table.querySelectorAll('th, td');
            cells.forEach(cell => {
                cell.style.padding = '12px 16px';
                cell.style.fontSize = '';
            });
        } else {
            table.classList.add('compact-mode');
            const cells = table.querySelectorAll('th, td');
            cells.forEach(cell => {
                cell.style.padding = '6px 8px';
                cell.style.fontSize = '0.9rem';
            });
        }

        // 重新调整列宽
        setTimeout(() => {
            adjustTableColumns(table);
        }, 100);
    }

    // ========== Add Table Interaction ==========
    function addTableInteraction(table, container) {

        // 行悬停效果
        const rows = table.querySelectorAll('tbody tr');
        rows.forEach(row => {
            row.addEventListener('mouseenter', function() {
                this.style.background = 'rgba(37, 99, 235, 0.05)';
                this.style.transform = 'translateX(2px)';
                this.style.boxShadow = '2px 0 8px rgba(37, 99, 235, 0.1)';
            });

            row.addEventListener('mouseleave', function() {
                this.style.background = '';
                this.style.transform = '';
                this.style.boxShadow = '';
            });
        });

        // 列悬停效果
        let currentColumn = -1;
        table.addEventListener('mouseover', function(e) {
            const cell = e.target.closest('th, td');
            if (!cell) return;

            const cellIndex = Array.from(cell.parentNode.children).indexOf(cell);

            if (cellIndex !== currentColumn) {
                // 清除之前的高亮
                clearColumnHighlight(table);

                // 高亮当前列
                highlightColumn(table, cellIndex);
                currentColumn = cellIndex;
            }
        });

        table.addEventListener('mouseleave', function() {
            clearColumnHighlight(table);
            currentColumn = -1;
        });

        // 双击列头自动调整该列宽度
        const headers = table.querySelectorAll('th');
        headers.forEach((header, index) => {
            header.addEventListener('dblclick', function() {
                autoFitSingleColumn(table, index);
                showToast(`🔧 Column ${index + 1} auto-fitted`);
            });

            // 添加工具提示
            header.title = 'Double-click to auto-fit this column';
            header.style.cursor = 'pointer';
        });
    }

    // ========== Column Highlighting ==========
    function highlightColumn(table, columnIndex) {
        const rows = table.querySelectorAll('tr');
        rows.forEach(row => {
            const cell = row.children[columnIndex];
            if (cell) {
                cell.style.background = 'rgba(37, 99, 235, 0.08)';
                cell.style.borderLeft = '2px solid #3b82f6';
                cell.style.borderRight = '2px solid #3b82f6';
            }
        });
    }

    function clearColumnHighlight(table) {
        const cells = table.querySelectorAll('th, td');
        cells.forEach(cell => {
            cell.style.background = cell.tagName === 'TH' ? 'linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%)' : '';
            cell.style.borderLeft = '1px solid #e2e8f0';
            cell.style.borderRight = '1px solid #e2e8f0';
        });
    }

    // ========== Auto Fit Single Column ==========
    function autoFitSingleColumn(table, columnIndex) {
        const rows = table.querySelectorAll('tr');
        let maxWidth = 0;

        rows.forEach(row => {
            const cell = row.children[columnIndex];
            if (cell) {
                // 临时移除宽度限制来测量内容
                const originalStyle = cell.style.cssText;
                cell.style.width = 'auto';
                cell.style.minWidth = 'auto';
                cell.style.maxWidth = 'none';

                const width = cell.scrollWidth;
                maxWidth = Math.max(maxWidth, width);

                // 恢复样式
                cell.style.cssText = originalStyle;
            }
        });

        // 更新 colgroup
        const colgroup = table.querySelector('colgroup');
        if (colgroup && colgroup.children[columnIndex]) {
            colgroup.children[columnIndex].style.width = (maxWidth + 20) + 'px';
        }

        // 添加动画效果
        rows.forEach(row => {
            const cell = row.children[columnIndex];
            if (cell) {
                cell.style.transform = 'scale(1.02)';
                setTimeout(() => {
                    cell.style.transform = '';
                }, 200);
            }
        });
    }

    // ========== Enhance Table Appearance ==========
    function enhanceTableAppearance(table) {
        // 添加表格包装器悬停效果
        const wrapper = table.closest('.responsive-table-wrapper');
        if (wrapper) {
            wrapper.addEventListener('mouseenter', function() {
                this.style.boxShadow = '0 8px 25px rgba(37, 99, 235, 0.2)';
                this.style.transform = 'translateY(-2px)';
            });

            wrapper.addEventListener('mouseleave', function() {
                this.style.boxShadow = '0 4px 15px rgba(37, 99, 235, 0.15)';
                this.style.transform = 'translateY(0)';
            });
        }

        // 添加斑马纹效果
        const rows = table.querySelectorAll('tbody tr');
        rows.forEach((row, index) => {
            if (index % 2 === 1) {
                row.style.backgroundColor = '#f8fafc';
            }
        });

        // 为长函数名添加特殊样式
        const cells = table.querySelectorAll('td, th');
        cells.forEach(cell => {
            const text = cell.textContent.trim();
            if (text.includes('skyborn.') && text.length > 25) {
                cell.style.fontFamily = 'Monaco, "Lucida Console", monospace';
                cell.style.fontSize = '0.9rem';
                cell.style.color = '#1e40af';
                cell.style.fontWeight = '500';

                // 添加函数名特殊标识
                if (cell.querySelector('a')) {
                    const link = cell.querySelector('a');
                    link.style.textDecoration = 'none';
                    link.style.borderBottom = '1px dotted #3b82f6';
                    link.style.transition = 'all 0.2s ease';

                    link.addEventListener('mouseenter', function() {
                        this.style.backgroundColor = 'rgba(59, 130, 246, 0.1)';
                        this.style.borderBottom = '1px solid #3b82f6';
                        this.style.transform = 'translateY(-1px)';
                    });

                    link.addEventListener('mouseleave', function() {
                        this.style.backgroundColor = '';
                        this.style.borderBottom = '1px dotted #3b82f6';
                        this.style.transform = '';
                    });
                }
            }
        });
    }

    // ========== Toast Notification System ==========
    function showToast(message, duration = 2000) {
        // 移除已存在的toast
        const existingToast = document.querySelector('.table-toast');
        if (existingToast) {
            existingToast.remove();
        }

        const toast = document.createElement('div');
        toast.className = 'table-toast';
        toast.textContent = message;
        toast.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 25%, #93c5fd 50%, #60a5fa 75%, #3b82f6 100%);
            color: white;
            padding: 0.75rem 1.5rem;
            border-radius: 12px;
            box-shadow: 0 8px 25px rgba(59, 130, 246, 0.3);
            z-index: 10000;
            font-size: 0.9rem;
            font-weight: 900;
            transform: translateX(100%);
            transition: transform 0.3s cubic-bezier(0.175, 0.885, 0.32, 1.275);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(59, 130, 246, 0.2);
        `;

        document.body.appendChild(toast);

        // 动画显示
        setTimeout(() => {
            toast.style.transform = 'translateX(0)';
        }, 10);

        // 自动隐藏
        setTimeout(() => {
            toast.style.transform = 'translateX(100%)';
            setTimeout(() => {
                if (toast.parentNode) {
                    toast.remove();
                }
            }, 300);
        }, duration);
    }

    // ========== Monitor DOM Changes ==========
    function monitorTableChanges() {
        const observer = new MutationObserver(function(mutations) {
            mutations.forEach(function(mutation) {
                if (mutation.type === 'childList') {
                    mutation.addedNodes.forEach(function(node) {
                        if (node.nodeType === Node.ELEMENT_NODE) {
                            // 检查新添加的表格
                            const newTables = node.querySelectorAll ? node.querySelectorAll('table') : [];
                            newTables.forEach(table => {
                                if (!table.classList.contains('responsive-table-processed')) {
                                    setTimeout(() => {
                                        makeTableResponsive(table);
                                    }, 100);
                                }
                            });

                            // 检查节点本身是否是表格
                            if (node.tagName === 'TABLE' && !node.classList.contains('responsive-table-processed')) {
                                setTimeout(() => {
                                    makeTableResponsive(node);
                                }, 100);
                            }
                        }
                    });
                }
            });
        });

        observer.observe(document.body, {
            childList: true,
            subtree: true
        });
    }

    // ========== Initialize System ==========
    function initializeResponsiveTablesSystem() {
        console.log('🚀 Initializing responsive tables system...');

        // 延迟初始化以确保所有内容加载完成
        setTimeout(() => {
            initializeTableAdjustment();
            monitorTableChanges();

            console.log('✅ Responsive tables system fully loaded!');

            // Create module-specific toast message
            const currentUrl = window.location.pathname;
            let moduleName = 'Smart tables';
            let icon = '📊';

            if (currentUrl.includes('functions_classes')) {
                if (window.location.hash.includes('calc')) { moduleName = 'calc'; icon = '🧮'; }
                else if (window.location.hash.includes('spharm')) { moduleName = 'spharm'; icon = '🔮'; }
                else if (window.location.hash.includes('gridfill')) { moduleName = 'gridfill'; icon = '🌐'; }
                else if (window.location.hash.includes('gradients')) { moduleName = 'gradients'; icon = '📈'; }
                else if (window.location.hash.includes('causality')) { moduleName = 'causality'; icon = '🔗'; }
                else if (window.location.hash.includes('windspharm')) { moduleName = 'windspharm'; icon = '�'; }

                if (moduleName !== 'Smart tables') {
                    moduleName = `${moduleName} functions`;
                }
            }

            showToast('🌟 Welcome to Skyborn!', 2000);
        }, 500);
    }

    // ========== Start the System ==========
    initializeResponsiveTablesSystem();

    // 为调试添加全局函数
    window.SkybornTables = {
        adjustAll: () => {
            document.querySelectorAll('table').forEach(adjustTableColumns);
            showToast('🔧 All tables adjusted');
        },
        resetAll: () => {
            document.querySelectorAll('table').forEach(resetTableColumns);
            showToast('🔄 All tables reset');
        },
        info: () => {
            const tables = document.querySelectorAll('table');
            console.log(`📊 Found ${tables.length} tables`);
            return tables;
        }
    };
});
