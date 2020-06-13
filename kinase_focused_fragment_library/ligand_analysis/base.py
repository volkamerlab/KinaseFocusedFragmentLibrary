"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This module contains basic classes for the analysis of the combinatorial library.
"""

from pathlib import Path


class ChemblPreparer:
    """
    Class to prepare the raw ChEMBL dataset for the analysis of the combinatorial library: Load raw data,
    """

    def __init__(self, path_chembl_raw, path_chembl_out):
        self.path_chembl_raw = Path(path_chembl_raw)
        self.path_chembl_out = Path(path_chembl_out)

    def prepare_chembl(in_file, out_file):
        """
        Prepare ChEMBL dataset, i.e. load, filter, standardize, convert data to InChI and  write (ChEMBL compound ID,
        standardized InChI) to output csv file.
        Write logging output to log file with same file name.

        Parameters
        ----------
        in_file : pathlib.Path
            Path to downloaded ChEMBL file.
        out_file : pathlib.Path
            Path to standardized output ChEMBL dataset file.
        """

        in_file = Path(in_file)
        out_file = Path(out_file)

        # read raw ChEMBL data
        molecules = _read_chembl_raw(in_file)

        # filter raw ChEMBL data
        molecules = _filter_chembl_raw(molecules)

        # standardize molecules
        molecules = _standardize(molecules)

        # get InChIs
        inchis = _get_inchis(molecules)
        print(inchis.head())

        # save data to file
        logger.info(f'Save to {out_file}...')
        inchis.to_csv(out_file)

    def _read_chembl_raw(path_to_chembl):
        """
        Read the raw ChEMBL data and drop keep only ChEMBL ID and canonical SMILES columns.

        Parameters
        ----------
        path_to_chembl : pathlib.Path
            Path to raw ChEMBL data.

        Returns
        -------
        pandas.Series
            SMILES (canonical_smiles) with ChEMBL ID as index (chembl_id).
        """

        # read data with columns: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key
        logger.info(f'Read {path_to_chembl}...')
        molecules = pd.read_csv(path_to_chembl, sep='\t')
        logger.info(f'Number of initial ChEMBL molecules: {molecules.shape[0]}')

        return molecules[['chembl_id', 'canonical_smiles']].set_index('chembl_id').squeeze()

    def _filter_chembl_raw(molecules):
        """
        Filter ChEMBL dataset: Extract largest molecule from mixtures, drop SMILES duplicates, add RDKit molecules column,
        and keep only molecules with more than 5 heavy atoms.

        Parameters
        ----------
        molecules : pandas.Series
            SMILES (canonical_smiles) with ChEMBL ID as index (chembl_id).

        Returns
        -------
        pandas.Series
            RDKit molecules (ROMol) with ChEMBL ID as index (chembl_id).
        """

        # Overwrite SMILES column with SMILES for longest molecule in original SMILES
        logger.info(f'Number of mixture SMILES: {molecules[molecules.str.contains(".", regex=False)].shape[0]}')
        molecules = molecules.apply(lambda x: max(x.split('.'), key=len))

        # Drop SMILES duplicates
        molecules.drop_duplicates(inplace=True)
        logger.info(f'Number of molecules after SMILES deduplication: {molecules.shape[0]}')

        # Add RDKit molecules column and drop SMILES column
        molecules = molecules.reset_index()
        PandasTools.AddMoleculeColumnToFrame(molecules, 'canonical_smiles')

        # Keep only molecules with more than 5 heavy atoms
        molecules['n_atoms'] = molecules.apply(lambda x: x.ROMol.GetNumHeavyAtoms(), axis=1)
        molecules = molecules[molecules.n_atoms > 5].copy()
        logger.info(f'Number of molecules after filter for number of atoms: {molecules.shape[0]}')

        return molecules[['chembl_id', 'ROMol']].set_index('chembl_id').squeeze()

    def _standardize(molecules):
        """
        Standardize molecules (ROMol).

        Parameters
        ----------
        molecules : pandas.Series
            RDKit molecules (ROMol) with ChEMBL ID as index (chembl_id).

        Returns
        -------
        pandas.Series
            Standardized RDKit molecules (ROMol) with ChEMBL ID as index (chembl_id).
        """

        # overwrite column with standardized molecules
        logger.info(f'Standardize molecules...')
        molecules = molecules.apply(standardize_mol)

        # drop rows with any data missing
        molecules.dropna(how='any', inplace=True)
        logger.info(f'Number of ChEMBL molecules after standardization: {molecules.shape[0]}')

        return molecules

    def _get_inchis(molecules):
        """

        Parameters
        ----------
        molecules : pandas.Series
            RDKit molecules (ROMol) with ChEMBL ID as index (chembl_id).

        Returns
        -------
        pandas.Series
            InChIs (inchi) with ChEMBL ID as index (chembl_id).
        """

        # molecules to InChIs
        inchis = molecules.apply(convert_mol_to_inchi)
        inchis.name = 'inchi'

        # drop rows with any data missing
        inchis.dropna(how='any', inplace=True)
        logger.info(f'Number of ChEMBL molecules (InChIs): {inchis.shape[0]}')

        return inchis


class CombinatorialLibraryAnalyzer():

    def __init__(self):
        pass