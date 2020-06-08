"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors
"""

# Add imports here
from . import preprocessing
from . import fragmentation
from . import recombination
from . import analysis

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions