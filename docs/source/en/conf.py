# English documentation configuration
# This configuration extends the main conf.py

import sys
import os
sys.path.insert(0, os.path.abspath('../../../src'))

# Import main configuration
exec(open('../conf.py').read())

# Language-specific overrides
language = 'en'
html_title = 'Skyborn Documentation'

# Master document for English
master_doc = 'index'

# Add any English-specific configuration here
