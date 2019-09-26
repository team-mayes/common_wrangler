"""
common_wrangler
Scripts for file and data manipulation, as well as common scripts used in this project as well as gaussian_wrangler and md_wrangler
"""

# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Add imports here
from .common import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
