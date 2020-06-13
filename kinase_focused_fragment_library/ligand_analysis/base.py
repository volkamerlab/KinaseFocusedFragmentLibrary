"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This module contains basic classes for the analysis of the combinatorial library.
"""


class ChemblPreparer:
    """
    Class to prepare the raw ChEMBL dataset for the analysis of the combinatorial library: Load raw data,
    """

    def __init__(self, path_chembl_raw, path_chembl_out):
        self.path_chembl_raw = path_chembl_raw
        self.path_chembl_out = path_chembl_out


class CombinatorialLibraryAnalyzer():

    def __init__(self):
        pass