/* Copy Button Fix for Skyborn Documentation */

// Wait for the DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {
    // Function to fix copy button positioning
    function fixCopyButtonPositioning() {
        // Find all copy buttons that might be misplaced
        const copyButtons = document.querySelectorAll('button.copybtn');

        copyButtons.forEach(button => {
            // Check if the button is not inside a highlight container
            const highlightContainer = button.closest('div.highlight');
            const parentHighlight = button.parentElement?.querySelector?.('div.highlight') ||
                                  button.nextElementSibling?.matches?.('div.highlight') ||
                                  button.previousElementSibling?.matches?.('div.highlight');

            if (!highlightContainer && parentHighlight) {
                // Move the button inside the highlight container
                parentHighlight.appendChild(button);
            } else if (!highlightContainer) {
                // Find the nearest highlight container
                const nearestHighlight = button.parentElement?.querySelector('div.highlight');
                if (nearestHighlight) {
                    nearestHighlight.appendChild(button);
                }
            }

            // Ensure proper positioning
            button.style.position = 'absolute';
            button.style.top = '0.3em';
            button.style.right = '0.3em';
            button.style.zIndex = '10';
        });
    }

    // Fix positioning immediately
    fixCopyButtonPositioning();

    // Fix positioning after a short delay (in case buttons are added dynamically)
    setTimeout(fixCopyButtonPositioning, 500);

    // Monitor for dynamically added content
    const observer = new MutationObserver(function(mutations) {
        let shouldFix = false;
        mutations.forEach(function(mutation) {
            if (mutation.type === 'childList') {
                mutation.addedNodes.forEach(function(node) {
                    if (node.nodeType === Node.ELEMENT_NODE) {
                        if (node.classList?.contains('copybtn') ||
                            node.querySelector?.('button.copybtn')) {
                            shouldFix = true;
                        }
                    }
                });
            }
        });

        if (shouldFix) {
            setTimeout(fixCopyButtonPositioning, 100);
        }
    });

    // Start observing
    observer.observe(document.body, {
        childList: true,
        subtree: true
    });
});

// Additional fallback: fix positioning on window load
window.addEventListener('load', function() {
    setTimeout(function() {
        const copyButtons = document.querySelectorAll('button.copybtn');
        copyButtons.forEach(button => {
            const highlightContainer = button.closest('div.highlight');
            if (highlightContainer) {
                // Ensure the highlight container has relative positioning
                highlightContainer.style.position = 'relative';
            }

            // Ensure button stays in top-right
            button.style.position = 'absolute';
            button.style.top = '0.3em';
            button.style.right = '0.3em';
            button.style.zIndex = '10';
        });
    }, 100);
});
