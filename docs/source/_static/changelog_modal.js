document.addEventListener('DOMContentLoaded', function () {
    if (typeof DOCUMENTATION_OPTIONS === 'undefined' ||
        DOCUMENTATION_OPTIONS.pagename !== 'changelog') {
        return;
    }

    const changelogRoot = document.querySelector('.bd-article > section#changelog');
    if (!changelogRoot || changelogRoot.dataset.modalReady === 'true') {
        return;
    }
    changelogRoot.dataset.modalReady = 'true';

    injectChangelogModalStyles();

    const versionSections = Array.from(changelogRoot.children).filter(function (child) {
        return child.tagName === 'SECTION' && child.id.indexOf('version-') === 0;
    });

    if (!versionSections.length) {
        return;
    }

    const versionMap = {};
    addChangelogHint(changelogRoot, versionSections[0]);

    versionSections.forEach(function (section) {
        const heading = getDirectHeading(section);
        if (!heading) {
            return;
        }

        const titleText = getHeadingText(heading);
        const statusText = getStatusText(titleText);
        const content = ensureSectionContent(section, heading);
        const modalButton = createModalButton(section, titleText);

        versionMap[section.id] = section;
        section.classList.add('changelog-modal-section');
        heading.classList.add('changelog-section-toggle');
        heading.tabIndex = 0;
        heading.setAttribute('role', 'button');
        heading.setAttribute('aria-label', 'Toggle release notes for ' + titleText);
        heading.appendChild(modalButton);
        content.classList.add('changelog-section-content');

        setSectionExpanded(section, isDefaultExpanded(statusText));

        heading.addEventListener('click', function (event) {
            if (event.target.closest('.changelog-modal-open-button') ||
                event.target.closest('.headerlink')) {
                return;
            }

            event.preventDefault();
            toggleSection(section);
        });

        heading.addEventListener('keydown', function (event) {
            if (event.target !== heading) {
                return;
            }

            if (event.key === 'Enter' || event.key === ' ') {
                event.preventDefault();
                toggleSection(section);
            }
        });
    });

    document.addEventListener('click', function (event) {
        if (!(event.target instanceof Element)) {
            return;
        }

        const link = event.target.closest('a[href^="#version-"]');
        if (!link) {
            return;
        }

        const targetId = (link.getAttribute('href') || '').slice(1);
        const targetSection = versionMap[targetId];
        if (!targetSection) {
            return;
        }

        event.preventDefault();
        expandSectionAndScroll(targetSection, targetId);
    });

    const hashTarget = window.location.hash ? window.location.hash.slice(1) : '';
    if (hashTarget && versionMap[hashTarget]) {
        setSectionExpanded(versionMap[hashTarget], true);
    }
});

function injectChangelogModalStyles() {
    if (document.getElementById('changelog-modal-styles')) {
        return;
    }

    const style = document.createElement('style');
    style.id = 'changelog-modal-styles';
    style.textContent = `
        .changelog-modal-tip {
            margin: 0 0 1.5rem;
            padding: 0.9rem 1.1rem;
            border-radius: 16px;
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.14), rgba(59, 130, 246, 0.08));
            border: 1px solid rgba(37, 99, 235, 0.18);
            color: #1e3a8a;
            font-weight: 600;
            box-shadow: 0 8px 24px rgba(37, 99, 235, 0.08);
        }

        .changelog-modal-section {
            margin: 0 0 1.5rem;
            padding: 0.45rem 1.1rem 1.15rem;
            border-radius: 18px;
            background: linear-gradient(180deg, rgba(255, 255, 255, 0.98), rgba(219, 234, 254, 0.58));
            border: 1px solid rgba(37, 99, 235, 0.12);
            box-shadow: 0 10px 28px rgba(37, 99, 235, 0.08);
        }

        .changelog-section-toggle {
            position: relative;
            display: block;
            cursor: pointer;
            margin: 0.25rem 0 0.85rem;
            padding: 0.95rem 9rem 0.95rem 3rem;
            border-radius: 16px;
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.14), rgba(59, 130, 246, 0.08));
            transition: transform 0.2s ease, box-shadow 0.2s ease, background 0.2s ease;
            outline: none;
        }

        .changelog-section-toggle::before {
            content: '';
            position: absolute;
            top: 50%;
            left: 1rem;
            width: 0;
            height: 0;
            border-top: 0.34rem solid transparent;
            border-bottom: 0.34rem solid transparent;
            border-left: 0.5rem solid #1d4ed8;
            transform: translateY(-50%) rotate(90deg);
            transform-origin: 35% 50%;
            transition: transform 0.2s ease;
        }

        .changelog-modal-section.is-collapsed .changelog-section-toggle::before {
            transform: translateY(-50%) rotate(0deg);
        }

        .changelog-section-toggle:hover,
        .changelog-section-toggle:focus {
            transform: translateY(-2px);
            box-shadow: 0 12px 28px rgba(37, 99, 235, 0.14);
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.2), rgba(59, 130, 246, 0.12));
        }

        .changelog-modal-open-button {
            position: absolute;
            top: 50%;
            right: 1rem;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-width: 7.25rem;
            transform: translateY(-50%);
            padding: 0.32rem 0.72rem;
            border: none;
            border-radius: 999px;
            background: rgba(37, 99, 235, 0.12);
            color: #1d4ed8;
            font-size: 0.78rem;
            font-weight: 700;
            letter-spacing: 0.02em;
            text-transform: uppercase;
            white-space: nowrap;
            cursor: pointer;
            transition: transform 0.18s ease, background 0.18s ease, color 0.18s ease;
        }

        .changelog-modal-open-button:hover,
        .changelog-modal-open-button:focus {
            transform: translateY(-50%) scale(1.03);
            background: rgba(37, 99, 235, 0.2);
            outline: none;
        }

        .changelog-section-content {
            padding: 0 0.3rem 0.1rem;
        }

        .changelog-modal-section.is-collapsed .changelog-section-content {
            display: none;
        }

        .changelog-modal-overlay {
            position: fixed;
            inset: 0;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 2rem;
            background: rgba(15, 23, 42, 0.72);
            backdrop-filter: blur(14px);
            z-index: 10000;
            opacity: 0;
            pointer-events: none;
            transition: opacity 0.22s ease;
        }

        .changelog-modal-overlay.is-visible {
            opacity: 1;
            pointer-events: auto;
        }

        .changelog-modal {
            width: min(980px, 96vw);
            max-height: 88vh;
            display: flex;
            flex-direction: column;
            background: linear-gradient(180deg, #ffffff 0%, #eff6ff 100%);
            border-radius: 24px;
            border: 1px solid rgba(255, 255, 255, 0.65);
            box-shadow: 0 24px 64px rgba(15, 23, 42, 0.28);
            overflow: hidden;
            transform: translateY(18px) scale(0.97);
            transition: transform 0.22s ease;
        }

        .changelog-modal-overlay.is-visible .changelog-modal {
            transform: translateY(0) scale(1);
        }

        .changelog-modal-open {
            overflow: hidden;
        }

        .changelog-modal-header {
            display: flex;
            align-items: flex-start;
            justify-content: space-between;
            gap: 1rem;
            padding: 1.4rem 1.5rem 1rem;
            border-bottom: 1px solid rgba(37, 99, 235, 0.12);
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.12), rgba(59, 130, 246, 0.06));
        }

        .changelog-modal-badge {
            display: inline-flex;
            align-items: center;
            margin: 0 0 0.55rem;
            padding: 0.28rem 0.65rem;
            border-radius: 999px;
            background: rgba(29, 78, 216, 0.12);
            color: #1d4ed8;
            font-size: 0.78rem;
            font-weight: 700;
            letter-spacing: 0.04em;
            text-transform: uppercase;
        }

        .changelog-modal-title {
            margin: 0;
            color: #0f172a;
            font-size: clamp(1.45rem, 2.5vw, 2rem);
            line-height: 1.15;
        }

        .changelog-modal-meta {
            margin: 0.5rem 0 0;
            color: #475569;
            font-size: 0.95rem;
        }

        .changelog-modal-topics {
            display: flex;
            flex-wrap: wrap;
            gap: 0.55rem;
            margin: 0 0 1.1rem;
        }

        .changelog-modal-topic {
            display: inline-flex;
            align-items: center;
            padding: 0.35rem 0.7rem;
            border-radius: 999px;
            background: rgba(37, 99, 235, 0.09);
            color: #1e3a8a;
            font-size: 0.86rem;
            font-weight: 600;
        }

        .changelog-modal-close {
            flex: 0 0 auto;
            width: 2.4rem;
            height: 2.4rem;
            border: none;
            border-radius: 999px;
            background: rgba(255, 255, 255, 0.88);
            color: #334155;
            font-size: 1rem;
            font-weight: 700;
            cursor: pointer;
            transition: transform 0.18s ease, background 0.18s ease, color 0.18s ease;
        }

        .changelog-modal-close:hover,
        .changelog-modal-close:focus {
            transform: scale(1.06);
            background: rgba(239, 68, 68, 0.12);
            color: #b91c1c;
            outline: none;
        }

        .changelog-modal-body {
            overflow: auto;
            padding: 1.35rem 1.5rem 1.5rem;
        }

        .changelog-modal-body p {
            color: #1e293b;
        }

        .changelog-modal-content > p {
            margin: 0 0 0.8rem;
            color: #1e3a8a;
            font-weight: 700;
            font-size: 1rem;
            letter-spacing: 0.01em;
        }

        .changelog-modal-content ul,
        .changelog-modal-content ol {
            margin: 0 0 1.15rem 1.25rem;
        }

        .changelog-modal-content li {
            margin-bottom: 0.5rem;
            color: #0f172a;
            line-height: 1.65;
        }

        [data-theme="dark"] .changelog-modal-tip {
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.24), rgba(59, 130, 246, 0.12));
            border-color: rgba(96, 165, 250, 0.2);
            color: #bfdbfe;
        }

        [data-theme="dark"] .changelog-modal-section {
            background: linear-gradient(180deg, rgba(15, 23, 42, 0.94), rgba(30, 41, 59, 0.88));
            border-color: rgba(96, 165, 250, 0.14);
            box-shadow: 0 14px 32px rgba(2, 6, 23, 0.35);
        }

        [data-theme="dark"] .changelog-section-toggle {
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.26), rgba(59, 130, 246, 0.14));
            color: #eff6ff !important;
        }

        [data-theme="dark"] .changelog-section-toggle::before {
            border-left-color: #bfdbfe;
        }

        [data-theme="dark"] .changelog-modal-open-button {
            background: rgba(96, 165, 250, 0.18);
            color: #bfdbfe;
        }

        [data-theme="dark"] .changelog-modal-open-button:hover,
        [data-theme="dark"] .changelog-modal-open-button:focus {
            background: rgba(96, 165, 250, 0.26);
        }

        [data-theme="dark"] .changelog-modal {
            background: linear-gradient(180deg, #0f172a 0%, #1e293b 100%);
            border-color: rgba(148, 163, 184, 0.16);
        }

        [data-theme="dark"] .changelog-modal-header {
            border-bottom-color: rgba(96, 165, 250, 0.16);
            background: linear-gradient(135deg, rgba(37, 99, 235, 0.22), rgba(59, 130, 246, 0.08));
        }

        [data-theme="dark"] .changelog-modal-title,
        [data-theme="dark"] .changelog-modal-content li,
        [data-theme="dark"] .changelog-modal-body p {
            color: #e2e8f0;
        }

        [data-theme="dark"] .changelog-modal-meta {
            color: #94a3b8;
        }

        [data-theme="dark"] .changelog-modal-content > p,
        [data-theme="dark"] .changelog-modal-topic {
            color: #bfdbfe;
        }

        [data-theme="dark"] .changelog-modal-topic {
            background: rgba(96, 165, 250, 0.12);
        }

        [data-theme="dark"] .changelog-modal-close {
            background: rgba(30, 41, 59, 0.88);
            color: #e2e8f0;
        }

        @media (max-width: 768px) {
            .changelog-modal-overlay {
                padding: 1rem;
            }

            .changelog-section-toggle {
                padding: 0.95rem 1rem 3.1rem 3rem;
            }

            .changelog-modal-open-button {
                top: auto;
                bottom: 0.85rem;
                right: 1rem;
                transform: none;
            }

            .changelog-modal-open-button:hover,
            .changelog-modal-open-button:focus {
                transform: scale(1.03);
            }

            .changelog-modal-header {
                padding: 1.1rem 1rem 0.85rem;
            }

            .changelog-modal-body {
                padding: 1rem;
            }
        }
    `;

    document.head.appendChild(style);
}

function addChangelogHint(changelogRoot, firstSection) {
    const hint = document.createElement('div');
    hint.className = 'changelog-modal-tip';
    hint.textContent = 'Planned and Current stay expanded by default. Click any version heading to expand or collapse it, and use Open modal for the popup view.';
    changelogRoot.insertBefore(hint, firstSection);
}

function getDirectHeading(section) {
    return Array.from(section.children).find(function (child) {
        return child.tagName === 'H2';
    }) || null;
}

function getSectionContent(section) {
    return section.querySelector(':scope > .changelog-section-content');
}

function ensureSectionContent(section, heading) {
    const existingContent = getSectionContent(section);
    if (existingContent) {
        return existingContent;
    }

    const content = document.createElement('div');
    content.className = 'changelog-section-content';

    Array.from(section.children).forEach(function (child) {
        if (child !== heading) {
            content.appendChild(child);
        }
    });

    section.appendChild(content);
    return content;
}

function getHeadingText(heading) {
    const clone = heading.cloneNode(true);
    clone.querySelectorAll('.headerlink, .changelog-modal-open-button').forEach(function (node) {
        node.remove();
    });
    return clone.textContent.replace(/\s+/g, ' ').trim();
}

function getStatusText(titleText) {
    const match = titleText.match(/\(([^)]+)\)\s*$/);
    return match ? match[1] : 'Release';
}

function isDefaultExpanded(statusText) {
    return statusText === 'Planned' || statusText === 'Current';
}

function createModalButton(section, titleText) {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'changelog-modal-open-button';
    button.textContent = 'Open modal';
    button.setAttribute('aria-label', 'Open modal for ' + titleText);

    button.addEventListener('click', function (event) {
        event.preventDefault();
        event.stopPropagation();
        openChangelogModal(section, button);
    });

    return button;
}

function setSectionExpanded(section, expanded) {
    const heading = getDirectHeading(section);
    const content = getSectionContent(section);

    section.classList.toggle('is-expanded', expanded);
    section.classList.toggle('is-collapsed', !expanded);

    if (heading) {
        heading.setAttribute('aria-expanded', expanded ? 'true' : 'false');
    }

    if (content) {
        content.hidden = !expanded;
    }
}

function toggleSection(section) {
    setSectionExpanded(section, !section.classList.contains('is-expanded'));
}

function expandSectionAndScroll(section, targetId) {
    setSectionExpanded(section, true);

    if (history.replaceState) {
        history.replaceState(null, '', '#' + targetId);
    } else {
        window.location.hash = targetId;
    }

    section.scrollIntoView({
        behavior: 'smooth',
        block: 'start'
    });
}

function getTopicLabels(section) {
    const content = getSectionContent(section);
    if (!content) {
        return [];
    }

    return Array.from(content.children)
        .filter(function (child) {
            return child.tagName === 'P' &&
                child.nextElementSibling &&
                (child.nextElementSibling.tagName === 'UL' ||
                    child.nextElementSibling.tagName === 'OL');
        })
        .map(function (child) {
            return child.textContent.trim();
        })
        .filter(Boolean)
        .slice(0, 6);
}

function openChangelogModal(section, opener) {
    const existingOverlay = document.querySelector('.changelog-modal-overlay');
    if (existingOverlay) {
        existingOverlay.remove();
    }

    const heading = getDirectHeading(section);
    const content = getSectionContent(section);
    if (!heading || !content) {
        return;
    }

    const titleText = getHeadingText(heading);
    const statusText = getStatusText(titleText);
    const detailCount = content.querySelectorAll('li').length;

    const overlay = document.createElement('div');
    overlay.className = 'changelog-modal-overlay';

    const modal = document.createElement('div');
    modal.className = 'changelog-modal';
    modal.setAttribute('role', 'dialog');
    modal.setAttribute('aria-modal', 'true');
    modal.setAttribute('aria-labelledby', 'changelog-modal-title');

    const header = document.createElement('div');
    header.className = 'changelog-modal-header';

    const headerText = document.createElement('div');

    const badge = document.createElement('span');
    badge.className = 'changelog-modal-badge';
    badge.textContent = statusText;

    const title = document.createElement('h2');
    title.id = 'changelog-modal-title';
    title.className = 'changelog-modal-title';
    title.textContent = titleText;

    const meta = document.createElement('p');
    meta.className = 'changelog-modal-meta';
    meta.textContent = detailCount + ' changelog item' + (detailCount === 1 ? '' : 's');

    headerText.appendChild(badge);
    headerText.appendChild(title);
    headerText.appendChild(meta);

    const closeButton = document.createElement('button');
    closeButton.type = 'button';
    closeButton.className = 'changelog-modal-close';
    closeButton.setAttribute('aria-label', 'Close release notes');
    closeButton.textContent = 'x';

    header.appendChild(headerText);
    header.appendChild(closeButton);

    const body = document.createElement('div');
    body.className = 'changelog-modal-body';

    const topics = getTopicLabels(section);
    if (topics.length) {
        const topicContainer = document.createElement('div');
        topicContainer.className = 'changelog-modal-topics';
        topics.forEach(function (label) {
            const tag = document.createElement('span');
            tag.className = 'changelog-modal-topic';
            tag.textContent = label;
            topicContainer.appendChild(tag);
        });
        body.appendChild(topicContainer);
    }

    const contentClone = document.createElement('div');
    contentClone.className = 'changelog-modal-content';

    Array.from(content.children).forEach(function (child) {
        contentClone.appendChild(child.cloneNode(true));
    });

    body.appendChild(contentClone);
    modal.appendChild(header);
    modal.appendChild(body);
    overlay.appendChild(modal);
    document.body.appendChild(overlay);
    document.body.classList.add('changelog-modal-open');

    requestAnimationFrame(function () {
        overlay.classList.add('is-visible');
    });

    function closeModal() {
        overlay.classList.remove('is-visible');
        document.body.classList.remove('changelog-modal-open');
        document.removeEventListener('keydown', escHandler);

        setTimeout(function () {
            if (overlay.parentNode) {
                overlay.remove();
            }
            if (opener && typeof opener.focus === 'function') {
                opener.focus();
            }
        }, 220);
    }

    function escHandler(event) {
        if (event.key === 'Escape') {
            closeModal();
        }
    }

    closeButton.addEventListener('click', closeModal);
    overlay.addEventListener('click', function (event) {
        if (event.target === overlay) {
            closeModal();
        }
    });
    document.addEventListener('keydown', escHandler);

    closeButton.focus();
}
