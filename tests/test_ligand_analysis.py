"""
Unit and regression test for the kinase_focused_fragment_library.ligand_analysis module.
"""

from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

from kinase_focused_fragment_library.ligand_analysis.base import ChemblPreparer

#PATH_TEST_DATA = Path(__name__).parent / 'kinase_focused_fragment_library' / 'tests' / 'data'  # TODO make work
PATH_TEST_DATA = Path('/home/dominique/Documents/Work/GitHub/KinaseFocusedFragmentLibrary/tests/data')


class TestsChemblPreparer:
    """
    Test ChemblPreparer class methods.
    """

    @pytest.mark.parametrize('path_chembl_raw, n_rows', [
        (
            PATH_TEST_DATA / 'chembl_25_chemreps_top100.txt',
            100
        )
    ])
    def test_read(self, path_chembl_raw, n_rows):
        """
        Test reading data.

        Parameters
        ----------
        path_chembl_raw : pathlib.Path or str
            Path to raw ChEMBL data (input).
        """

        chembl_preparer = ChemblPreparer()
        molecules = chembl_preparer._read(path_chembl_raw=path_chembl_raw)

        assert isinstance(molecules, pd.Series)
        assert molecules.name == 'canonical_smiles'
        assert molecules.index.name == 'chembl_id'
        assert molecules.shape[0] == n_rows

    @pytest.mark.parametrize('smiles, n_rows', [
        (
            pd.DataFrame(
                {
                    'chembl_id': [
                        'CHEMBL153534',
                        'CHEMBL501667'
                    ],
                    'canonical_smiles': [
                        'Cc1cc(cn1C)c2csc(N=C(N)N)n2',
                        '[Br-].[Br-].CCn1c2ccc3cc2c4cc(ccc14)C(=O)c5ccc(Cn6cc[n+](Cc7ccc(cc7)c8cccc(c9ccc(C[n+]%10ccn(Cc%11ccc(cc%11)C3=O)c%10)cc9)c8C(=O)O)c6)cc5'
                    ]
                }
            ).set_index('chembl_id').squeeze(),
            2
        )
    ])
    def test_filter(self, smiles, n_rows):
        """
        Test filtering data.

        Parameters
        ----------
        smiles : pandas.Series
            SMILES ("canonical_smiles") with ChEMBL ID as index ("chembl_id").
        n_rows : int
            Number of DataFrame rows.
        """

        chembl_preparer = ChemblPreparer()
        molecules_filtered = chembl_preparer._filter(smiles=smiles)

        assert isinstance(molecules_filtered, pd.Series)
        assert molecules_filtered.name == 'ROMol'
        assert molecules_filtered.index.name == 'chembl_id'
        assert molecules_filtered.shape[0] == n_rows

    @pytest.mark.parametrize('molecules, n_rows', [
        (
            pd.DataFrame(
                {
                    'chembl_id': [
                        'CHEMBL153534',
                        'CHEMBL501667'
                    ],
                    'canonical_smiles': [
                        'Cc1cc(cn1C)c2csc(N=C(N)N)n2',
                        'CCn1c2ccc3cc2c4cc(ccc14)C(=O)c5ccc(Cn6cc[n+](Cc7ccc(cc7)c8cccc(c9ccc(C[n+]%10ccn(Cc%11ccc(cc%11)C3=O)c%10)cc9)c8C(=O)O)c6)cc5'
                    ]
                }
            ).set_index('chembl_id').squeeze(),
            2
        )
    ])
    def test_standardize(self, molecules, n_rows):
        """
        Test standardizing data.

        Parameters
        ----------
        molecules : pandas.Series
            RDKit molecules ("ROMol") with ChEMBL ID as index ("chembl_id").
        n_rows : int
            Number of DataFrame rows.
        """

        chembl_preparer = ChemblPreparer()

        # convert SMILES to molecules
        molecules = molecules.apply(lambda x: Chem.MolFromSmiles(x))
        molecules.name = 'ROMol'

        molecules_standardized = chembl_preparer._standardize(molecules=molecules)

        assert isinstance(molecules_standardized, pd.Series)
        assert molecules_standardized.name == 'ROMol'
        assert molecules_standardized.index.name == 'chembl_id'
        assert molecules_standardized.shape[0] == n_rows

    @pytest.mark.parametrize('molecules, n_rows', [
        (
            pd.DataFrame(
                {
                    'chembl_id': [
                        'CHEMBL153534',
                        'CHEMBL501667'
                    ],
                    'canonical_smiles': [
                        'Cc1cc(cn1C)c2csc(N=C(N)N)n2',
                        'CCn1c2ccc3cc2c4cc(ccc14)C(=O)c5ccc(Cn6cc[n+](Cc7ccc(cc7)c8cccc(c9ccc(C[n+]%10ccn(Cc%11ccc(cc%11)C3=O)c%10)cc9)c8C(=O)O)c6)cc5'
                    ]
                }
            ).set_index('chembl_id').squeeze(),
            2
        )
    ])
    def test_get_inchis(self, molecules, n_rows):
        """
        Test getting InChIs.

        Parameters
        ----------
        molecules : pandas.Series
            RDKit molecules ("ROMol") with ChEMBL ID as index ("chembl_id").
        n_rows : int
            Number of DataFrame rows.
        """

        chembl_preparer = ChemblPreparer()

        # convert SMILES to molecules
        molecules = molecules.apply(lambda x: Chem.MolFromSmiles(x))
        molecules.name = 'ROMol'

        inchis = chembl_preparer._get_inchis(molecules=molecules)

        assert isinstance(inchis, pd.Series)
        assert inchis.name == 'inchi'
        assert inchis.index.name == 'chembl_id'
        assert inchis.shape[0] == n_rows
