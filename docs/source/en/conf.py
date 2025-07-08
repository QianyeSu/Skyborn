# -*- coding: utf-8 -*-
# English documentation configuration
# This configuration extends the main conf.py

import sys
import os
sys.path.insert(0, os.path.abspath('../../../src'))

# Import main configuration with proper encoding
with open('../conf.py', 'r', encoding='utf-8') as f:
    exec(f.read())

# Language-specific overrides
language = 'en'
html_title = 'Skyborn Documentation'

# Master document for English
master_doc = 'index'

# Add any English-specific configuration here
