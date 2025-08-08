/* Skyborn Documentation - Dynamic Table Column Width Adjustment */

document.addEventListener('DOMContentLoaded', function() {

    // ========== Table Auto-Adjustment System ==========
    function initializeTableAdjustment() {

        // Find all tables
        const tables = document.querySelectorAll('table');

        tables.forEach(table => {
            makeTableResponsive(table);
        });

        // Listen for window size changes - but don't auto-adjust manually adjusted tables
        let resizeTimer;
        window.addEventListener('resize', function() {
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(() => {
                tables.forEach(table => {
                    // Only adjust tables that haven't been manually adjusted by user
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

        // Skip already processed tables
        if (table.classList.contains('responsive-table-processed')) {
            return;
        }

        // Add responsive container
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

        // Create table container
        const container = document.createElement('div');
        container.className = 'table-scroll-container';
        container.style.cssText = `
            width: 100%;
            overflow-x: auto;
            overflow-y: hidden;
            position: relative;
            scroll-behavior: smooth;
            /* Custom scrollbar styles */
            scrollbar-width: thin;
            scrollbar-color: #3b82f6 #f1f5f9;
        `;

        // WebKit browser scrollbar styles
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

        // Wrap table
        table.parentNode.insertBefore(wrapper, table);
        wrapper.appendChild(container);
        container.appendChild(table);

        // Add table controls
        addTableControls(wrapper, table, container);

        // Initial table adjustment - but respect RST defined column widths
        applyInitialLayout(table);

        // Add interaction effects
        addTableInteraction(table, container);

        // Mark as processed
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

        // Left side: Table title (obtained from preceding heading)
        const info = document.createElement('div');
        const tableTitleData = getTableTitle(table);

        // Create clickable title
        const titleElement = document.createElement('div');
        titleElement.innerHTML = tableTitleData.title;
        titleElement.style.cssText = `
            font-weight: 700;
            cursor: pointer;
            color: #1e40af;
            font-size: ${getHeadingSize(tableTitleData.headingLevel)};
            margin: 0;
            padding: 0.25rem 0;
            transition: all 0.3s ease;
            border-bottom: 2px solid transparent;
            user-select: none;
            position: relative;
        `;

        // Add hover indicator (# symbol)
        const hoverIndicator = document.createElement('span');
        hoverIndicator.innerHTML = ' #';
        hoverIndicator.style.cssText = `
            color: #3b82f6;
            font-size: 0.9em;
            opacity: 0;
            transition: all 0.3s ease;
            margin-left: 0.3rem;
        `;
        titleElement.appendChild(hoverIndicator);

        // Click to jump to Quick Navigation
        titleElement.addEventListener('click', () => {
            // Find Quick Navigation element
            const quickNavigation = document.getElementById('quick-navigation');
            if (quickNavigation) {
            quickNavigation.scrollIntoView({ behavior: 'smooth', block: 'start' });
            // Highlight Quick Navigation
            quickNavigation.style.transition = 'background-color 0.5s ease';
            quickNavigation.style.backgroundColor = 'rgba(59, 130, 246, 0.1)';
            quickNavigation.style.borderRadius = '8px';
            quickNavigation.style.padding = '1rem';
            setTimeout(() => {
                quickNavigation.style.backgroundColor = '';
                quickNavigation.style.padding = '';
            }, 2000);
            // Show success notification
            showToast('🧭 Returned to Quick Navigation');
            } else {
            // If Quick Navigation not found, jump to page top
            document.querySelector('h1').scrollIntoView({ behavior: 'smooth', block: 'start' });
            showToast('Returned to page top');
            }
        });

        titleElement.addEventListener('mouseenter', () => {
            titleElement.style.color = '#3b82f6';
            titleElement.style.borderBottomColor = '#3b82f6';
            titleElement.style.transform = 'scale(1.02)';
            titleElement.title = 'Click to return to Quick Navigation';
            hoverIndicator.style.opacity = '1';
            hoverIndicator.style.transform = 'translateY(-1px)';
        });

        titleElement.addEventListener('mouseleave', () => {
            titleElement.style.color = '#1e40af';
            titleElement.style.borderBottomColor = 'transparent';
            titleElement.style.transform = 'scale(1)';
            hoverIndicator.style.opacity = '0';
            hoverIndicator.style.transform = 'translateY(0)';
        });

        info.appendChild(titleElement);

        // Hide original title
        if (tableTitleData.headingElement) {
            tableTitleData.headingElement.style.display = 'none';
        }

        controlBar.appendChild(info);

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
        // Find table heading in front
        let currentElement = table.parentElement;
        let targetHeading = null;

        // Search up for heading in container
        while (currentElement && currentElement !== document.body) {
            // Find preceding heading elements
            const prevElements = [];
            let prev = currentElement.previousElementSibling;

            // Collect several preceding elements
            while (prev && prevElements.length < 5) {
                prevElements.unshift(prev);
                prev = prev.previousElementSibling;
            }

            // Search for heading in preceding elements
            for (const element of prevElements) {
                // Check if it's a heading element
                if (element.tagName && /^H[1-6]$/.test(element.tagName)) {
                    targetHeading = element;
                    return {
                        title: element.textContent.trim(),
                        headingElement: element,
                        headingLevel: element.tagName
                    };
                }

                // Check if contains heading text
                if (element.classList.contains('section-header') ||
                    element.classList.contains('table-title')) {
                    return {
                        title: element.textContent.trim(),
                        headingElement: element,
                        headingLevel: 'H3'
                    };
                }

                // Find heading in child elements
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

        // If no heading found, infer from table content
        const firstHeaderCell = table.querySelector('th');
        if (firstHeaderCell) {
            const headerText = firstHeaderCell.textContent.trim();

            // Infer table type based on first column heading
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

        // Check if table is in specific module section
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

        // Default title
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

        // Apply preset 30:70 ratio for functions_classes tables
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

            // Apply to all rows
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
            // Apply slight auto-adjustment for other tables
            adjustTableColumns(table);
        }
    }

    // ========== Adjust Table Columns ==========
    function adjustTableColumns(table) {
        const rows = table.querySelectorAll('tr');
        if (rows.length === 0) return;

        // Get first row to determine column count
        const firstRow = rows[0];
        const cells = firstRow.querySelectorAll('th, td');
        const columnCount = cells.length;

        if (columnCount === 0) return;

        // Calculate maximum content width for each column
        const columnWidths = new Array(columnCount).fill(0);
        const columnContents = new Array(columnCount).fill().map(() => []);

        rows.forEach(row => {
            const cells = row.querySelectorAll('th, td');
            cells.forEach((cell, index) => {
                if (index < columnCount) {
                    const textLength = cell.textContent.trim().length;
                    const hasLongFunction = cell.textContent.includes('skyborn.') && textLength > 25;

                    // Allocate more space for long function names
                    const adjustedLength = hasLongFunction ? textLength * 1.2 : textLength;
                    columnWidths[index] = Math.max(columnWidths[index], adjustedLength);
                    columnContents[index].push(cell.textContent.trim());
                }
            });
        });

        // Calculate available container width
        const container = table.closest('.table-scroll-container');
        const availableWidth = container ? container.clientWidth - 40 : window.innerWidth - 100;

        // Smart column width allocation
        const totalContentWidth = columnWidths.reduce((sum, width) => sum + width, 0);
        const avgCharWidth = 8; // Estimated character width
        const estimatedTableWidth = totalContentWidth * avgCharWidth;

        // If content width exceeds available space, use proportional allocation
        if (estimatedTableWidth > availableWidth) {
            const scale = availableWidth / estimatedTableWidth;
            columnWidths.forEach((width, index) => {
                columnWidths[index] = Math.max(width * scale, 100); // Minimum width 100px
            });
        }

        // Apply column widths
        applyColumnWidths(table, columnWidths);

        // Add table appearance enhancement
        enhanceTableAppearance(table);
    }

    // ========== Apply Column Widths ==========
    function applyColumnWidths(table, widths) {
        // Set table styles
        table.style.cssText = `
            width: 100%;
            table-layout: fixed;
            border-collapse: collapse;
            margin: 0;
            background: white;
            transition: all 0.3s ease;
        `;

        // Create or update colgroup
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

        // Apply cell styles
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

        // Header styles
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
        // Temporarily remove table-layout: fixed to measure content
        table.style.tableLayout = 'auto';
        table.style.width = 'auto';

        // Force browser to recalculate
        table.offsetHeight;

        setTimeout(() => {
            // Get natural widths
            const rows = table.querySelectorAll('tr');
            const firstRow = rows[0];
            const cells = firstRow.querySelectorAll('th, td');
            const naturalWidths = [];

            cells.forEach(cell => {
                naturalWidths.push(cell.offsetWidth);
            });

            // Reapply fixed layout
            applyColumnWidths(table, naturalWidths);

            // Add animation effect
            table.style.transform = 'scale(1.01)';
            setTimeout(() => {
                table.style.transform = 'scale(1)';
            }, 200);

        }, 50);
    }

    // ========== Reset Table Columns ==========
    function resetTableColumns(table) {
        // Remove custom styles
        table.style.cssText = '';

        const colgroup = table.querySelector('colgroup');
        if (colgroup) {
            colgroup.remove();
        }

        // Reinitialize
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

        // Readjust column widths
        setTimeout(() => {
            adjustTableColumns(table);
        }, 100);
    }

    // ========== Add Table Interaction ==========
    function addTableInteraction(table, container) {

        // Row hover effects
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

        // Column hover effects
        let currentColumn = -1;
        table.addEventListener('mouseover', function(e) {
            const cell = e.target.closest('th, td');
            if (!cell) return;

            const cellIndex = Array.from(cell.parentNode.children).indexOf(cell);

            if (cellIndex !== currentColumn) {
                // Clear previous highlight
                clearColumnHighlight(table);

                // Highlight current column
                highlightColumn(table, cellIndex);
                currentColumn = cellIndex;
            }
        });

        table.addEventListener('mouseleave', function() {
            clearColumnHighlight(table);
            currentColumn = -1;
        });

        // Double-click column header to auto-adjust that column width
        const headers = table.querySelectorAll('th');
        headers.forEach((header, index) => {
            header.addEventListener('dblclick', function() {
                autoFitSingleColumn(table, index);
                showToast(`🔧 Column ${index + 1} auto-fitted`);
            });

            // Add tooltip
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
                // Temporarily remove width restrictions to measure content
                const originalStyle = cell.style.cssText;
                cell.style.width = 'auto';
                cell.style.minWidth = 'auto';
                cell.style.maxWidth = 'none';

                const width = cell.scrollWidth;
                maxWidth = Math.max(maxWidth, width);

                // Restore styles
                cell.style.cssText = originalStyle;
            }
        });

        // Update colgroup
        const colgroup = table.querySelector('colgroup');
        if (colgroup && colgroup.children[columnIndex]) {
            colgroup.children[columnIndex].style.width = (maxWidth + 20) + 'px';
        }

        // Add animation effect
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
        // Add table wrapper hover effects
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

        // Add zebra stripe effects
        const rows = table.querySelectorAll('tbody tr');
        rows.forEach((row, index) => {
            if (index % 2 === 1) {
                row.style.backgroundColor = '#f8fafc';
            }
        });

        // Add special styles for long function names
        const cells = table.querySelectorAll('td, th');
        cells.forEach(cell => {
            const text = cell.textContent.trim();
            if (text.includes('skyborn.') && text.length > 25) {
                cell.style.fontFamily = 'Monaco, "Lucida Console", monospace';
                cell.style.fontSize = '0.9rem';
                cell.style.color = '#1e40af';
                cell.style.fontWeight = '500';

                // Add special identifier for function names
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
        // Remove existing toast
        const existingToast = document.querySelector('.table-toast');
        if (existingToast) {
            existingToast.remove();
        }

        // Detect dark mode
        const isDarkMode = document.documentElement.getAttribute('data-theme') === 'dark' ||
                          window.matchMedia('(prefers-color-scheme: dark)').matches ||
                          document.body.classList.contains('dark') ||
                          document.body.classList.contains('dark-mode');

        const toast = document.createElement('div');
        toast.className = 'table-toast';
        toast.textContent = message;

        // Choose gradient based on mode
        const lightGradient = 'linear-gradient(135deg, #e0f2fe 0%, #7dd3fc 25%, #38bdf8 50%, #0ea5e9 75%, #0284c7 100%)';
        const darkGradient = 'linear-gradient(135deg, #1e293b 0%, #334155 25%, #475569 50%, #1e40af 75%, #1d4ed8 100%)';

        toast.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            background: ${isDarkMode ? darkGradient : lightGradient};
            color: ${isDarkMode ? '#e2e8f0' : 'white'};
            padding: 0.75rem 1.5rem;
            border-radius: 12px;
            box-shadow: 0 8px 25px rgba(59, 130, 246, 0.3);
            z-index: 10000;
            font-size: 0.9rem;
            font-weight: 750;
            transform: translateX(100%);
            transition: transform 0.3s cubic-bezier(0.175, 0.885, 0.32, 1.275);
            backdrop-filter: blur(10px);
            border: 1px solid ${isDarkMode ? 'rgba(100, 116, 139, 0.3)' : 'rgba(59, 130, 246, 0.2)'};
        `;

        document.body.appendChild(toast);

        // Animation show
        setTimeout(() => {
            toast.style.transform = 'translateX(0)';
        }, 10);

        // Auto hide
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
                            // Check new added tables
                            const newTables = node.querySelectorAll ? node.querySelectorAll('table') : [];
                            newTables.forEach(table => {
                                if (!table.classList.contains('responsive-table-processed')) {
                                    setTimeout(() => {
                                        makeTableResponsive(table);
                                    }, 100);
                                }
                            });

                            // Check if node itself is a table
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

        // Delayed initialization to ensure all content is loaded
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
                else if (window.location.hash.includes('windspharm')) { moduleName = 'windspharm'; icon = '💨'; }

                if (moduleName !== 'Smart tables') {
                    moduleName = `${moduleName} functions`;
                }
            }

            showToast('🌟 Welcome to Skyborn!', 3000);
        }, 500);
    }

    // ========== Start the System ==========
    initializeResponsiveTablesSystem();

    // Add global functions for debugging
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
