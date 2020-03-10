"""
Utility functions to work with the fragment library.
"""

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina


def read_fragment_library(path_to_lib, remove_dummy=True, reduced=''):
    """
    Read fragment library from sdf files (one file per subpocket).

    Parameters
    ----------
    path_to_lib : str
        Path to fragment library folder.
    remove_dummy : bool
        Replace dummy atoms with hydrogens in fragments (default), or leave dummy atoms in fragments.
    reduced : str
        Add here string in file name describing the library reduction, e.g. add '_reduced_0.5' if file name is
        'AP_reduced_0.5.sdf'.


    Returns
    -------
    dict of pandas.DataFrame
        Fragment details, i.e. SMILES, and fragment RDKit molecules, KLIFS and fragmentation details (values)
        for each subpocket (key).
    """

    # list of folders for each subpocket
    subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2', 'X']

    data = {}

    # iterate over subpockets
    for subpocket in subpockets:

        if subpocket=='X':
            data[subpocket] = _read_subpocket_fragments(subpocket, path_to_lib, remove_dummy)
        else:
            data[subpocket] = _read_subpocket_fragments(subpocket, path_to_lib, remove_dummy, reduced)

    return data


def _read_subpocket_fragments(subpocket, path_to_lib, remove_dummy=True, reduced=''):
    """
    Read fragments for input subpocket.

    Parameters
    ----------
    subpocket : str
        Subpocket name, i.e. AP, SE, FP, GA, B1, or B2.
    path_to_lib : str
        Path to fragment library folder.
    remove_dummy : bool
        Replace dummy atoms with hydrogens in fragments (default), or leave dummy atoms in fragments.
    reduced : str
        Add here string in file name describing the library reduction, e.g. add '_reduced_0.5' if file name is
        'AP_reduced_0.5.sdf'.

    Returns
    -------
    pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, for input subpocket.
    """

    try:
        mol_supplier = Chem.SDMolSupplier(str(path_to_lib / f'{subpocket}{reduced}.sdf'), removeHs=False)
    except OSError:
        mol_supplier = Chem.SDMolSupplier(str(path_to_lib / subpocket / f'{subpocket}{reduced}.sdf'), removeHs=False)

    data = {
        'smiles': [],
        'fragment': [],
        'kinase': [],
        'family': [],
        'group': [],
        'complex_pdb': [],
        'ligand_pdb': [],
        'alt': [],
        'chain': [],
        'klifs_code': [],
        'atom.prop.subpocket': [],
        'atom.prop.environment': []
    }

    for mol_raw in mol_supplier:

        if remove_dummy:

            # Replace dummy atoms with hydrogens in fragments
            dummy = Chem.MolFromSmiles('*')
            hydrogen = Chem.MolFromSmiles('[H]', sanitize=False)
            mol = AllChem.ReplaceSubstructs(mol_raw, dummy, hydrogen, replaceAll=True)[0]
            mol = Chem.RemoveHs(mol)  # Remove all hydrogens but explicit hydrogens

        else:

            mol = Chem.RemoveHs(mol_raw)

        # Generate SMILES for comparing fragments
        smiles = Chem.MolToSmiles(mol)
        data['smiles'].append(smiles)

        # 2D coordinates
        AllChem.Compute2DCoords(mol)
        data['fragment'].append(mol)

        # KLIFS data
        data['kinase'].append(mol.GetProp('kinase'))
        data['family'].append(mol.GetProp('family'))
        data['group'].append(mol.GetProp('group'))
        data['complex_pdb'].append(mol.GetProp('complex_pdb'))
        data['ligand_pdb'].append(mol.GetProp('ligand_pdb'))
        data['alt'].append(mol.GetProp('alt'))
        data['chain'].append(mol.GetProp('chain'))
        data['klifs_code'].append(mol.GetProp('_Name'))

        # Fragmentation data
        data['atom.prop.subpocket'].append(mol.GetProp('atom.prop.subpocket'))
        data['atom.prop.environment'].append(mol.GetProp('atom.prop.environment'))

    return pd.DataFrame(data)


def get_deduplicated_original_complexes(fragment_library):
    """
    Get the set of complexes (kinase name and KLIFS complex code) from which the fragment library originates.

    Parameters
    ----------
    fragment_library : dict of pandas.DataFrame
        Fragment details, i.e. SMILES, and fragment RDKit molecules, KLIFS and fragmentation details (values)
        for each subpocket (key).

    Returns
    -------
    pandas.DataFrame
        Deduplicated kinase-structure pairs.
    """

    # Concatenate all subpockets to one DataFrame
    fragment_library_concat = pd.concat([values for key, values in fragment_library.items()], axis=0)
    print(f'Fragments: {len(fragment_library_concat)}')

    # Get deduplicated kinase-structure pairs
    kinase_structure_pairs = pd.DataFrame(
        {
            'kinase': fragment_library_concat.kinase,
            'klifs_code': fragment_library_concat.klifs_code
        }
    )
    kinase_structure_pairs.drop_duplicates(inplace=True)
    print(f'Deduplicated kinase-structure pairs from which fragments originate: '
          f'{len(kinase_structure_pairs)}')

    return kinase_structure_pairs


def most_common_fragments(fragments, top_x=50):
    """
    Get most common fragments.

    Parameters
    ----------
    fragments : pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, for input subpocket.

    top_x : int
        Top x most common fragments.

    Returns
    -------
    tuple (list of rdkit.Chem.rdchem.Mol, pandas.Series)
        List of top x fragments (RDKit molecules) and frequence of top x fragments in subpocket (Series).
    """

    # Sort fragments by number of counts
    mols_count = fragments.smiles.value_counts()  # Sorted in descending order

    # Get RDKit Mol from SMILES
    mols = [Chem.MolFromSmiles(smiles) for smiles in mols_count.index]

    # N most common fragments
    return mols[:top_x], mols_count[:top_x]


def generate_fingerprints(mols):
    """
    Generate RDKit fingerprint from list of molecules.

    Parameters
    ----------
    mols : list of rdkit.Chem.rdchem.Mol
        List of molecules.

    Returns
    -------
    list of rdkit.DataStructs.cDataStructs.ExplicitBitVect
        List of fingerprints.
    """

    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    fingerprints = [rdkit_gen.GetFingerprint(mol) for mol in mols]

    return fingerprints


def cluster_molecules(fingerprints, cutoff=0.6):
    """
    Cluster fingerprints using the Butina algorithm.

    Parameters
    ----------
    fingerprints : list of rdkit.DataStructs.cDataStructs.ExplicitBitVect
        List of fingerprints.
    cutoff : float
        Distance cutoff Butina clustering.

    Returns
    -------
    list of tuple of int
        List of clusters, whereby each cluster is described by its cluster member IDs.
    """

    # Calculate Tanimoto distance matrix
    distance_matrix = _tanimoto_distance_matrix(fingerprints)

    # Now cluster the data with the implemented Butina algorithm:
    clusters = Butina.ClusterData(
        distance_matrix,
        len(fingerprints),
        cutoff,
        isDistData=True
    )

    # Sort clusters by size
    clusters = sorted(clusters, key=len, reverse=True)

    # Get number of singleton clusters
    num_singletons = len([cluster for cluster in clusters if len(cluster) == 1])

    # Print details on clustering
    print("Number of fragments:", len(fingerprints))
    print("Threshold: ", cutoff)
    print("Number of clusters: ", len(clusters))
    print("# clusters with only 1 compound: ", num_singletons)

    return clusters


def _tanimoto_distance_matrix(fingerprints):
    """
    Calculate distance matrix for list of fingerprints.

    Parameters
    ----------
    fingerprints : list of rdkit.DataStructs.cDataStructs.ExplicitBitVect
        List of fingerprints.

    Returns
    -------
    list of floats
        Distance matrix (a triangular similarity matrix in the form of a list)
    """

    fingerprints = list(fingerprints)
    distance_matrix = []

    for i in range(1, len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[:i])

        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        distance_matrix.extend([1 - x for x in similarities])

    return distance_matrix
