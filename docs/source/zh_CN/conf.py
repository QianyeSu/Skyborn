# Chinese documentation configuration
# This configuration extends the main conf.py

import sys
import os
sys.path.insert(0, os.path.abspath('../../../src'))

# Import main configuration
exec(open('../conf.py').read())

# Language-specific overrides
language = 'zh_CN'
html_title = 'Skyborn 文档'

# Master document for Chinese
master_doc = 'index'

# Add any Chinese-specific configuration here
