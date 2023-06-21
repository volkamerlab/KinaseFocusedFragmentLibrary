"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors
"""

# Add imports here
import kinase_focused_fragment_library.preprocessing
import kinase_focused_fragment_library.fragmentation
import kinase_focused_fragment_library.fragment_library_reduction
import kinase_focused_fragment_library.recombination
import kinase_focused_fragment_library.ligand_analysis

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
