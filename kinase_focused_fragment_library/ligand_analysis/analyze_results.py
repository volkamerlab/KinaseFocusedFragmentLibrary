import json
import multiprocessing as mp
import time

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

from .utils import standardize_mol, construct_ligand
from kinase_focused_fragment_library.recombination.pickle_loader import pickle_loader


def get_ligands_analysis(fragment_library, original_ligands, chembl, in_paths, combinatorial_library_file):
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
    in_paths : list of pathlib.Path
        Paths to ligand pickle files.
    combinatorial_library_file : pathlib.Path
        Path to output json file containing combinatorial library meta data and properties.
    """

    start = time.time()

    n_processes = 2  # mp.cpu_count()
    print("Number of processors: ", n_processes)
    pool = mp.Pool(n_processes)

    results = []

    # iterate over pickle files
    for in_path in in_paths:

        print(str(in_path))
        with open(str(in_path), 'rb') as pickle_in:

            # process ligands in pickle file (returns list of dict)
            results_tmp = pool.starmap(
                _get_ligand_analysis,
                [(meta, fragment_library, original_ligands, chembl) for meta in pickle_loader(pickle_in)]
            )
            print(f'Number of ligands in current iteration: {len(results_tmp)}')

            # extend results list with ligands from current iteration
            results.extend(results_tmp)

    print(f'Number of ligands from all iterations: {len(results)}')

    with open(str(combinatorial_library_file), 'w') as f:
        json.dump(results, f)

    runtime = time.time() - start
    print('Time: ', runtime)


def _get_ligand_analysis(meta, fragment_library, original_ligands, chembl):
    """
    Construct ligand from fragment library based on meta data (fragment and bond ids) and analyze ligand with respect
    to the following properties:
    - Lipinski's rule of five properties: HBA, HBD, molecular weight, and logP
    - Number of atoms
    - Exact matches in ChEMBL molecule dataset
    - Exact matches or substructure matches in original ligand dataset

    Parameters
    ----------
    meta : TBA
        TBA
    fragment_library : dict of pandas.DataFrame
        Fragment library, i.e. fragments (value) per subpocket (key).
    original_ligands : pandas.DataFrame
        Standardized original ligands (ligands from with fragment library is originating): InCHI and ROMol.
    chembl : pandas.DataFrame
        Standardized ChEMBL molecules: InCHI.

    Returns
    -------
    dict
        Ligand's fragment IDs, bond IDs, Lipinski's rule of five criteria, exact matches in ChEMBL and in the original
        ligand as well as substructure matches in the original ligands.
    """

    print('hallo!')

    # initialize ligand details
    ligand_dict = {
        'bond_ids': [list(i) for i in meta.bonds],
        'fragment_ids': list(meta.frag_ids),
        'hba': None,
        'hbd': None,
        'mwt': None,
        'logp': None,
        'n_atoms': None,
        'chembl_exact': None,
        'original_exact': None,
        'original_substructure': None,
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
    except Exception as e:
        print(f'Error {e}: Mol to InChI conversion failed for molecule with fragment IDs {result.meta.frag_ids}')
        return

    # Lipinski's rule of five
    lipinski, wt, logp, hbd, hba = is_drug_like(ligand)

    # number of atoms
    n_atoms = ligand.GetNumHeavyAtoms()

    # ligand has exact match in original ligands?
    original_exact_matches = original_ligands[
        original_ligands.inchi == inchi
    ].index.to_list()

    # ligand has substructure match in original ligands?
    original_substructure_matches = original_ligands[
        original_ligands.mol.apply(lambda x: x.HasSubstructMatch(ligand))
    ].index.to_list()

    # ligand has exact match in ChEMBL?
    chembl_exact_matches = chembl[
        chembl.standard_inchi == inchi
    ].index.to_list()

    # save results to dictionary
    ligand_dict['hba'] = hba
    ligand_dict['hbd'] = hbd
    ligand_dict['mwt'] = wt
    ligand_dict['logp'] = logp
    ligand_dict['n_atoms'] = n_atoms
    ligand_dict['chembl_exact'] = chembl_exact_matches
    ligand_dict['original_exact'] = original_exact_matches
    ligand_dict['original_substructure'] = original_substructure_matches
    ligand_dict['inchi'] = inchi

    return ligand_dict


def is_drug_like(mol):
    """
    Get Lipinski's rule of five criteria for molecule.

    (If used in loop for multiple molecules, it takes about 1s for 2000 molecules.)

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    tuple of int
        Fulfilled criteria (1) or not (0) for Lipinski's rule of five, and its criteria molecule weight, logP, HBD and
        HBA.
    """

    mol_wt = 1 if Descriptors.ExactMolWt(mol) <= 500 else 0

    logp = 1 if Descriptors.MolLogP(mol) <= 5 else 0

    hbd = 1 if Lipinski.NumHDonors(mol) <= 5 else 0

    hba = 1 if Lipinski.NumHAcceptors(mol) <= 10 else 0

    lipinski = 1 if mol_wt + logp + hbd + hba >= 3 else 0

    return lipinski, mol_wt, logp, hbd, hba