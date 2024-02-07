"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This module contains basic classes for the analysis of the combinatorial library.
"""

import json
import multiprocessing as mp
import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import PandasTools, rdFingerprintGenerator

from .utils import standardize_smiles_to_inchi, standardize_mol, construct_ligand, is_drug_like
from kinase_focused_fragment_library.recombination.pickle_loader import pickle_loader

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

        # filter raw ChEMBL data
        smiles = self._filter(smiles)

        # standardize InChIs
        inchis = self._standardize(smiles)

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

        logger.info(f'Read raw ChEMBL data...')
        logger.info(f'Path: {path_chembl_raw}')

        # read data with columns: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key
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
            RDKit molecules ("canonical_smiles") with ChEMBL ID as index ("chembl_id").
        """

        logger.info(f'Filter molecules (SMILES)...')

        # get SMILES for longest molecule in original SMILES
        logger.info(f'Number of mixture SMILES: {smiles[smiles.str.contains(".", regex=False)].shape[0]}')
        smiles = smiles.apply(lambda x: max(x.split('.'), key=len))

        # drop SMILES duplicates
        smiles.drop_duplicates(inplace=True)
        logger.info(f'Number of molecules after SMILES deduplication: {smiles.shape[0]}')

        # get RDKit molecules (conversion returns None if molecules cannot be constructed)
        molecules = smiles.reset_index()
        PandasTools.AddMoleculeColumnToFrame(molecules, 'canonical_smiles')

        # drop rows with any data missing
        molecules.dropna(how='any', inplace=True)
        logger.info(f'Number of ChEMBL molecules after dropping None values (ROMol): {molecules.shape[0]}')

        # keep only molecules with more than 5 heavy atoms
        molecules['n_atoms'] = molecules.apply(lambda x: x.ROMol.GetNumHeavyAtoms(), axis=1)
        molecules = molecules[molecules.n_atoms > 5]
        logger.info(f'Number of molecules after filter for number of atoms: {molecules.shape[0]}')

        # cast to Series
        smiles_filtered = molecules[['chembl_id', 'canonical_smiles']].set_index('chembl_id').squeeze()

        return smiles_filtered

    @staticmethod
    def _standardize(smiles):
        """
        Standardize molecules (from SMILES to InChi).

        Parameters
        ----------
        smiles : pandas.Series
            Molecule SMILES ("canonical_smiles") with ChEMBL ID as index ("chembl_id").

        Returns
        -------
        pandas.Series
            Standardized molecule InChIs ("inchi") with ChEMBL ID as index ("chembl_id").
        """

        logger.info(f'Standardize molecules (ROMol)...')

        # get standardized molecules
        inchis = smiles.apply(standardize_smiles_to_inchi)
        inchis.name = 'inchi'

        # drop rows with any data missing
        inchis.dropna(how='any', inplace=True)
        logger.info(f'Number of ChEMBL InChIs after standardization: {inchis.shape[0]}')

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

        logger.info(f'Save molecules (InChI)...')
        logger.info(f'Path: {path_chembl_out}')
        inchis.to_csv(path_chembl_out)


class CombinatorialLibraryAnalyzer:

    def __init__(self):
        pass

    def run(self, fragment_library, original_ligands, chembl, path_combinatorial_library):
        """
        Construct ligands from fragment library based on meta data (fragment and bond ids) and analyze ligands with respect
        to the following properties:
        - Lipinski's rule of five properties: HBA, HBD, molecular weight, and logP
        - Number of atoms
        - Exact matches in ChEMBL molecule dataset
        - Exact matches or substructure matches in original ligand dataset

        Parameters
        ----------
        fragment_library : dict of pandas.DataFrame
            Fragment library, i.e. fragments (value) per subpocket (key).
        original_ligands : pandas.DataFrame
            Standardized original ligands (ligands from with fragment library is originating): InCHI and ROMol.
        chembl : pandas.DataFrame
            Standardized ChEMBL molecules: InCHI.
        path_combinatorial_library : pathlib.Path
            Path to combinatorial library folder, contains pickled molecules from recombination and is used to output final
            json file.
        """

        n_cores = mp.cpu_count()
        logger.info(f'Number of cores: {n_cores}')
        pool = mp.Pool(n_cores)

        results = []

        # get all paths to pickle files with molecules from recombination
        paths_pickle_combinatorial_library = list((path_combinatorial_library / 'results').glob('*.pickle'))

        logger.info(f'Process recombined ligands...')

        # iterate over pickle files
        for i, path_pickle_combinatorial_library in enumerate(paths_pickle_combinatorial_library):
            logger.info(f'Process {path_pickle_combinatorial_library}...')

            with open(str(path_pickle_combinatorial_library), 'rb') as pickle_in:
                # process ligands in pickle file (returns list of dict)
                results_tmp = pool.starmap(
                    self._analyze_ligand,
                    [(meta, fragment_library, original_ligands, chembl) for meta in pickle_loader(pickle_in)]
                )
                logger.info(f'Number of ligands in current iteration: {len(results_tmp)}')

                # extend results list with ligands from current iteration
                results.extend(results_tmp)

            with open(path_pickle_combinatorial_library.parent / f'{path_pickle_combinatorial_library.stem}.json', 'w') as f:
                json.dump(results_tmp, f)

        logger.info(f'Number of recombined ligands from all iterations: {len(results)}')
        logger.info(f'Data linked to each ligand: {list(results[0].keys())}')

        with open(path_combinatorial_library / 'combinatorial_library.json', 'w') as f:
            json.dump(results, f)

    def _analyze_ligand(self, meta, fragment_library, original_ligands, chembl):
        """
        Construct ligand from fragment library based on meta data (fragment and bond ids) and analyze ligand with respect
        to the following properties:
        - Lipinski's rule of five properties: HBA, HBD, molecular weight, and logP
        - Number of atoms
        - Exact matches in ChEMBL molecule dataset
        - Exact matches or substructure matches in original ligand dataset

        Parameters
        ----------
        meta : kinase_focused_fragment_library.recombination.classes_meta.Combination
            Ligand's meta data: fragment IDs (list of str) and bond IDs (list of list of str), where the strings are
            composed of: "subpocket_name"_"fragment_index".
        fragment_library : dict of pandas.DataFrame
            Fragment library, i.e. fragments (value) per subpocket (key).
        original_ligands : pandas.DataFrame
            Standardized original ligands (ligands from with fragment library is originating): InCHI and ROMol.
        chembl : pandas.DataFrame
            Standardized ChEMBL molecules: InCHI.

        Returns
        -------
        dict
            Ligand's fragment IDs, bond IDs, Lipinski's rule of five criteria, exact matches in ChEMBL, the highest
            Tanimoto similarity value to ChEMBL, and exact/substructure matches in the original ligands.
        """

        # initialize ligand details
        ligand_dict = {
            'bond_ids': [list(i) for i in meta.bonds],
            'fragment_ids': list(meta.frag_ids),
            'hba': None,
            'hbd': None,
            'mwt': None,
            'logp': None,
            'n_atoms': None,
            'chembl_exact': 0,
            'chembl_most_similar': None,
            'original_exact': 0,
            'original_substructure': 0,
            'inchi': None
        }

        # construct ligand
        try:
            ligand = construct_ligand(meta, fragment_library)
        except Exception as e:
            print(f'Error {e}: Molecule construction failed for molecule with fragment IDs {meta.frag_ids}')
            return

        # standardize molecule
        try:
            ligand = standardize_mol(ligand)
        except Exception as e:
            print(f'Error {e}: Molecule standardization failed for molecule with fragment IDs {meta.frag_ids}')
            return

        # convert mol to inchi
        try:
            inchi = Chem.MolToInchi(ligand)
            ligand_dict['inchi'] = inchi
        except Exception as e:
            print(f'Error {e}: Mol to InChI conversion failed for molecule with fragment IDs {meta.frag_ids}')
            return

        # Lipinski's rule of five
        lipinski, wt, logp, hbd, hba = is_drug_like(ligand)
        ligand_dict['hba'] = hba
        ligand_dict['hbd'] = hbd
        ligand_dict['mwt'] = wt
        ligand_dict['logp'] = logp

        # number of atoms
        ligand_dict['n_atoms'] = ligand.GetNumHeavyAtoms()

        # ligand has exact match in original ligands?
        if not original_ligands[original_ligands.inchi == inchi].empty:
            ligand_dict['original_exact'] = 1

        # ligand has substructure match in original ligands?
        if not original_ligands[original_ligands.mol.apply(lambda x: x.HasSubstructMatch(ligand))].empty:
            ligand_dict['original_substructure'] = 1

        # ligand has exact match in ChEMBL?
        if not chembl[chembl.inchi == inchi].empty:
            ligand_dict['chembl_exact'] = 1

        # highest Tanimoto similarity between ligand and ChEMBL?
        chembl_most_similar = self._most_similar_chembl_ligand(inchi, chembl)
        ligand_dict['chembl_most_similar'] = chembl_most_similar

        return ligand_dict

    @staticmethod
    def _most_similar_chembl_ligand(ligand_inchi, chembl):
        """
        Get the most similar ChEMBL ligand (ChEMBL compound ID and Tanimoto similarity) to the query ligand.

        Parameters
        ----------
        ligand_inchi : str
            Recombined ligand (InChI)
        chembl : pandas.DataFrame
            ChEMBL ligands, column fingerprint necessary.

        Returns
        -------
        tuple of (str, float)
            ChEMBL compound ID and Tanimoto similarity of ChEMBL ligand most similar to the query ligand.
        """

        try:

            # get ROMol from recombined ligand InChI
            ligand = Chem.MolFromInchi(ligand_inchi)

            # generate query ligand fingerprint
            rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
            query_fingerprint = rdkit_gen.GetFingerprint(ligand)

            # get ChEMBL fingerprints as list
            chembl_fingerprints = chembl.fingerprint.to_list()

            # get pairwise similarities
            chembl['similarity'] = DataStructs.BulkTanimotoSimilarity(query_fingerprint, chembl_fingerprints)

            # get ligand with maximal similarity
            chembl_most_similar_ix = chembl.similarity.idxmax()

            return [
                chembl.loc[chembl_most_similar_ix].chembl_id,
                round(chembl.loc[chembl_most_similar_ix].similarity, 2)
            ]

        except Exception as e:
            
            print(f'Most similar ChEMBL ligand search problem for {ligand_inchi}: {e}')
            return [None, None]


class CombinatorialLibraryDeduplicator:

    def __init__(self):
        pass

    @staticmethod
    def run(path_combinatorial_library):
        """
        Deduplicate combinatorial library based on InChIs.

        Parameters
        ----------
        path_combinatorial_library : pathlib.Path
            Path to combinatorial library folder, containing final combinatorial library json file
        """

        print('Load json file...')
        combinatorial_library = pd.read_json(path_combinatorial_library / 'combinatorial_library.json')
        print(f'Number of recombined ligands before deduplication: {combinatorial_library.shape[0]}')
        combinatorial_library.drop_duplicates('inchi', inplace=True)
        print(f'Number of recombined ligands after deduplication: {combinatorial_library.shape[0]}')
        print('Convert DataFrame to dict...')
        combinatorial_library = combinatorial_library.to_dict('records')
        print('Save as json file...')
        with open(path_combinatorial_library / 'combinatorial_library_deduplicated.json', 'w') as f:
            json.dump(combinatorial_library, f)
