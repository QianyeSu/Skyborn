# -- Path setup --
import datetime
import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Skyborn'
copyright = '2025, Qianye Su'
author = 'Qianye Su'

# Add last updated date
html_last_updated_fmt = '%Y-%m-%d'

# Dynamic version import with fallback
try:
    # Try to import skyborn to get version
    import skyborn
    release = skyborn.__version__
except (ImportError, AttributeError) as e:
    # Fallback: read version from pyproject.toml or __init__.py
    try:
        import toml
        with open('../../pyproject.toml', 'r') as f:
            pyproject = toml.load(f)
            release = pyproject['project']['version']
    except:
        # Final fallback: read from __init__.py
        import re
        with open('../../src/skyborn/__init__.py', 'r') as f:
            content = f.read()
            match = re.search(
                r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', content)
            release = match.group(1) if match else '0.3.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',  # 添加 autosummary 扩展
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx_autodoc_typehints',
    'myst_nb',  # For Jupyter notebook support (includes myst_parser)
    'sphinx_copybutton',  # Add copy button to code blocks
]

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = True

# Support for multiple file types (moved to end of file)

templates_path = ['_templates']
exclude_patterns = []

# Suppress specific warnings
suppress_warnings = [
    'myst.mathjax',
    # 'autodoc.duplicate_object',
    # 'autodoc',  # Add more general autodoc suppression
]

# -- Internationalization ---------------------------------------------------
language = 'en'  # Default language
locale_dirs = ['locale/']   # Path to locale directoryfuro
gettext_compact = False     # Generate POT files with long names

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

# Custom CSS and JS files
html_css_files = [
    'custom.css',
]

html_js_files = [
    'interactive.js',
    'copybutton_fix.js',
    'latex_formula_manager.js',
    'table_responsive.js',
]

# Theme options for sphinx_book_theme
html_theme_options = {
    "repository_url": "https://github.com/QianyeSu/skyborn",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_source_button": True,
    "use_issues_button": True,
    "use_download_button": True,
    "use_fullscreen_button": True,
    "path_to_docs": "docs/source",
    "repository_branch": "main",
    "launch_buttons": {
        "notebook_interface": "jupyterlab",
        "binderhub_url": "https://mybinder.org"
    },
    "show_navbar_depth": 2,
    "show_toc_level": 2,
    "navigation_with_keys": False,
    "collapse_navigation": False,
    "use_sidenotes": True,
    "home_page_in_toc": False,
}

# Set the logo link manually
html_logo_link = "index.html"

# -- Extension configuration -------------------------------------------------
# Autosummary configuration
autosummary_generate = True
autosummary_imported_members = True
autosummary_ignore_module_all = False

# Autodoc configuration
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
}

# -- Mock imports to avoid dependency issues --------------------------------
autodoc_mock_imports = [
    'metpy',
    'metpy.calc',
    'scipy',
    'scipy.ndimage',
    'scipy.special',
    'sklearn',
    'sklearn.feature_selection',
    'eccodes',
    'cfgrib',
    'netCDF4',
    'cartopy',
    'matplotlib',
    'seaborn',
    'xarray',
    'skyborn.ROF',  # ROF module is under development
]

# -- MyST-NB configuration for Jupyter notebooks ---------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# MyST-NB specific settings
nb_execution_mode = "off"  # Don't execute notebooks during build
nb_execution_timeout = 300  # Timeout for notebook execution (if enabled)
nb_execution_allow_errors = True  # Allow notebooks with errors
nb_merge_streams = True  # Merge stdout and stderr

# Support for both .md and .ipynb files
source_suffix = {
    '.rst': None,
    '.md': 'myst-nb',
    '.ipynb': 'myst-nb',
}

# MyST configuration to avoid conflicts
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    "substitution",
]

# -- MathJax configuration for better LaTeX support ------------------------
mathjax3_config = {
    "tex": {
        "inlineMath": [['$', '$'], ['\\(', '\\)']],
        "displayMath": [['$$', '$$'], ['\\[', '\\]']],
        "processEscapes": True,
        "processEnvironments": True,
        "macros": {
            "bm": [r"\boldsymbol{#1}", 1],
            "vec": [r"\mathbf{#1}", 1],
            "mat": [r"\mathbf{#1}", 1],
            "d": r"\mathrm{d}",
            "e": r"\mathrm{e}",
            "i": r"\mathrm{i}",
            "R": r"\mathbb{R}",
            "C": r"\mathbb{C}",
            "N": r"\mathbb{N}",
            "Z": r"\mathbb{Z}",
            "Q": r"\mathbb{Q}",
        },
        "packages": ["base", "ams", "newcommand", "configmacros"],
    },
    "options": {
        "ignoreHtmlClass": "tex2jax_ignore",
        "processHtmlClass": "tex2jax_process",
    },
    "loader": {
        "load": ["[tex]/ams", "[tex]/newcommand", "[tex]/configmacros"]
    },
    "svg": {
        "fontCache": "global"
    }
}

# Configure notebook output and MIME priorities
nb_output_stderr = "show"
nb_mime_priority_overrides = [
    ("html", "application/vnd.jupyter.widget-view+json", 10),
    ("html", "application/javascript", 20),
    ("html", "text/html", 30),
    ("html", "image/svg+xml", 40),
    ("html", "image/png", 50),
    ("html", "image/jpeg", 60),
    ("html", "text/markdown", 70),
    ("html", "text/latex", 80),
    ("html", "text/plain", 90),
]

# Improve notebook rendering
nb_render_image_options = {
    "align": "center",
    "width": "100%",
    "class": "responsive-image"
}

nb_render_figure_options = {
    "align": "center",
    "width": "100%",
    "class": "responsive-figure"
}

# -- HTML configuration for notebooks ---------------------------------------
html_title = f"{project} v{release}"
html_logo = "_static/SkyBornLogo.svg"
html_favicon = None

# Show source links for notebooks
html_show_sourcelink = True
html_copy_source = True

# Enable last updated date display
html_show_sphinx = True
html_show_copyright = True

# Configure copy button
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
