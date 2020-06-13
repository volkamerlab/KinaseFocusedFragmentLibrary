"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This module contains basic classes for the analysis of the combinatorial library.
"""

import logging
from pathlib import Path

import pandas as pd
from rdkit.Chem import PandasTools

from .utils import standardize_mol, convert_mol_to_inchi

logger = logging.getLogger(__name__)


class ChemblPreparer:
    """
    Class to prepare the raw ChEMBL dataset for the analysis of the combinatorial library.
    """

    def __init__(self):

        pass

    def run(self, path_chembl_raw, path_chembl_out):
        """
        Prepare ChEMBL dataset: Load, filter, standardize, and convert data to InChI
        Write ChEMBL compound ID and standardized InChI to output csv file, alongside a logging file.

        Parameters
        ----------
        path_chembl_raw : pathlib.Path or str
            Path to raw ChEMBL data (input).
        path_chembl_out : pathlib.Path or str
            Path to standardized ChEMBL data (output).
        """

        # set paths
        path_chembl_raw = Path(path_chembl_raw)
        path_chembl_out = Path(path_chembl_out)

        # read raw ChEMBL data
        smiles = self._read(path_chembl_raw)
        smiles = smiles[:1000]

        # filter raw ChEMBL data
        molecules = self._filter(smiles)

        # standardize molecules
        molecules = self._standardize(molecules)

        # get InChIs
        inchis = self._get_inchis(molecules)

        # save data to file
        self._save(inchis, path_chembl_out)

    @staticmethod
    def _read(path_chembl_raw):
        """
        Read the raw ChEMBL data and keep only ChEMBL ID and canonical SMILES.

        Parameters
        ----------
        path_chembl_raw : pathlib.Path or str
            Path to raw ChEMBL data (input).

        Returns
        -------
        pandas.Series
            SMILES ("canonical_smiles") with ChEMBL ID as index ("chembl_id").
        """

        # read data with columns: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key
        logger.info(f'Read {path_chembl_raw}...')
        molecules = pd.read_csv(path_chembl_raw, sep='\t')
        logger.info(f'Number of initial ChEMBL molecules: {molecules.shape[0]}')

        return molecules[['chembl_id', 'canonical_smiles']].set_index('chembl_id').squeeze()

    @staticmethod
    def _filter(smiles):
        """
        Filter ChEMBL dataset: Extract largest molecule from mixtures (canonical SMILES), drop SMILES duplicates,
        convert SMILES to RDKit molecules, and keep only molecules with more than 5 heavy atoms.

        Parameters
        ----------
        smiles : pandas.Series
            SMILES ("canonical_smiles") with ChEMBL ID as index ("chembl_id").

        Returns
        -------
        pandas.Series
            RDKit molecules ("ROMol") with ChEMBL ID as index ("chembl_id").
        """

        # get SMILES for longest molecule in original SMILES
        logger.info(f'Number of mixture SMILES: {smiles[smiles.str.contains(".", regex=False)].shape[0]}')
        smiles = smiles.apply(lambda x: max(x.split('.'), key=len))

        # drop SMILES duplicates
        smiles.drop_duplicates(inplace=True)
        logger.info(f'Number of molecules after SMILES deduplication: {smiles.shape[0]}')

        # get RDKit molecules
        molecules = smiles.reset_index()
        PandasTools.AddMoleculeColumnToFrame(molecules, 'canonical_smiles')

        # keep only molecules with more than 5 heavy atoms
        molecules['n_atoms'] = molecules.apply(lambda x: x.ROMol.GetNumHeavyAtoms(), axis=1)
        molecules = molecules[molecules.n_atoms > 5].copy()
        logger.info(f'Number of molecules after filter for number of atoms: {molecules.shape[0]}')

        return molecules[['chembl_id', 'ROMol']].set_index('chembl_id').squeeze()

    @staticmethod
    def _standardize(molecules):
        """
        Standardize molecules (ROMol).

        Parameters
        ----------
        molecules : pandas.Series
            RDKit molecules ("ROMol") with ChEMBL ID as index ("chembl_id").

        Returns
        -------
        pandas.Series
            Standardized RDKit molecules ("ROMol") with ChEMBL ID as index ("chembl_id").
        """

        # get standardized molecules
        logger.info(f'Standardize molecules...')
        molecules = molecules.apply(standardize_mol)

        # drop rows with any data missing
        molecules.dropna(how='any', inplace=True)
        logger.info(f'Number of ChEMBL molecules after standardization: {molecules.shape[0]}')

        return molecules

    @staticmethod
    def _get_inchis(molecules):
        """
        Convert RDKit molecules ("ROMol") to InChIs ("inchi").

        Parameters
        ----------
        molecules : pandas.Series
            RDKit molecules ("ROMol") with ChEMBL ID as index ("chembl_id").

        Returns
        -------
        pandas.Series
            InChIs ("inchi") with ChEMBL ID as index ("chembl_id").
        """

        # convert molecules to InChIs
        inchis = molecules.apply(convert_mol_to_inchi)
        inchis.name = 'inchi'

        # drop rows with any data missing
        inchis.dropna(how='any', inplace=True)
        logger.info(f'Number of ChEMBL molecules (InChIs): {inchis.shape[0]}')

        return inchis

    @staticmethod
    def _save(inchis, path_chembl_out):
        """
        Save ChEMBL IDs and InChIs to csv file.

        Parameters
        ----------
        inchis : pandas.Series
            InChIs ("inchi") with ChEMBL ID as index ("chembl_id").
        path_chembl_out : pathlib.Path or str
            Path to standardized ChEMBL data (output).
        """

        logger.info(f'Save to {path_chembl_out}...')
        inchis.to_csv(path_chembl_out)


class CombinatorialLibraryAnalyzer:

    def __init__(self):
        pass
