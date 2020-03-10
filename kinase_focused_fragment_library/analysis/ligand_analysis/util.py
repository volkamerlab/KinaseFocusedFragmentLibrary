"""
Utility functions to work with the combinatorial library.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def load_combinatorial_library_properties(recombined_ligand_properties_path, ligand_indices_deduplicated=None):
    """
    Load properties of combinatorial library (full or deduplicated entries).

    Parameters
    ----------
    recombined_ligand_properties_path : pathlib.Path
        Path to file with recombined ligand properties.
    ligand_indices_deduplicated : None or numpy.ndarray.
        Deduplicated ligand indices. Default is None, i.e. entries will not be deduplicated.

    Returns
    -------
    pandas.DataFrame
        Recombined ligand properties (full or deduplicated entries).
    """

    properties = np.genfromtxt(recombined_ligand_properties_path, delimiter=',', dtype=int)
    properties = pd.DataFrame(
        properties, columns=[
            'ligand_id',
            'mw',
            'hba',
            'hbd',
            'logp',
            'lipinski',
            'n_atoms',
            'n_subpockets',
            'original',
            'original_sub',
            'chembl_match'
        ]
    )
    properties.set_index('ligand_id', inplace=True)

    print(f'Full library: {properties.shape[0]}')

    if ligand_indices_deduplicated is not None:

        # Load ligand indices in deduplicated ligand library
        print(f'Ligand indices of deduplicated library entries: {len(ligand_indices_deduplicated)}')

        # Get only properties for deduplicated ligands
        properties = properties.loc[list(ligand_indices_deduplicated), :]
        print(f'Deduplicated library: {properties.shape[0]}')

    return properties


def load_ligand_indices_deduplicated(ligand_indices_deduplicated_path):
    """
    Load ligand indices of deduplicated combinatorial library from file.

    Parameters
    ----------
    ligand_indices_deduplicated_path : pathlib.Path
        Path to file with ligand indices.

    Returns
    -------
    numpy.ndarray
        Ligand indices of deduplicated combinatorial library.
    """

    ligand_indices_deduplicated_path = Path(ligand_indices_deduplicated_path)

    data = np.genfromtxt(ligand_indices_deduplicated_path, delimiter='/n', dtype=int)

    return data


def _drug_likeness_from_mol(mol):
    """
    Get drug-likeness criteria for a molecule, i.e. molecular weight, logP, number of hydrogen bond acceptors/donors and
    accordance to Lipinski's rule of five.
    (Takes about 1s for 2000 mols.)

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    pd.Series
        Drug-likeness criteria for input molecule.
    """

    mw = 1 if Descriptors.ExactMolWt(mol) <= 500 else 0
    logp = 1 if Descriptors.MolLogP(mol) <= 5 else 0
    hbd = 1 if Lipinski.NumHDonors(mol) <= 5 else 0
    hba = 1 if Lipinski.NumHAcceptors(mol) <= 10 else 0
    lipinski = 1 if mw + logp + hbd + hba >= 3 else 0

    return pd.Series([mw, logp, hbd, hba, lipinski], index='mw logp hbd hba lipinski'.split())


def drug_likeness_from_smiles(smiles):
    """
    Get drug-likeness for a set of SMILES.

    Parameters
    ----------
    smiles : pd.Series
        Set of SMILES.

    Returns
    -------
    pd.Series
        Ratio of drug like molecules.
    """

    drug_likeness = pd.DataFrame(
        smiles.apply(
            lambda x: _drug_likeness_from_mol(Chem.MolFromSmiles(x))
        )
    )
    print(f'Number of molecules: {drug_likeness.shape[0]}')

    drug_likeness_ratio = round(drug_likeness.apply(sum) / len(drug_likeness) * 100, 0)

    return drug_likeness_ratio


def drug_likeness_from_mols(mols):
    """
    Get drug-likeness for a set of molecule objects.

    Parameters
    ----------
    mols : list of rdkit.Chem.rdchem.Mol
        Molecules.

    Returns
    -------
    pd.Series
        Ratio of drug like molecules.
    """

    drug_likeness = pd.DataFrame(
        [_drug_likeness_from_mol(mol) for mol in mols]
    )
    print(f'Number of molecules: {drug_likeness.shape[0]}')

    drug_likeness_ratio = round(drug_likeness.apply(sum) / len(drug_likeness) * 100, 0)

    return drug_likeness_ratio


def drug_likeness_from_combinatorial_library(recombined_ligands_properties):
    """
    Get drug-likeness for combinatorial library.

    Parameters
    ----------
    recombined_ligands_properties : pandas.DataFrame
        Recombined ligand properties (full or deduplicated entries).

    Returns
    -------
    pd.Series
        Ratio of drug like molecules.
    """

    # Number of ligands that fullfill the druglikeness criteria
    drug_likeness = recombined_ligands_properties[['mw', 'logp', 'hbd', 'hba', 'lipinski']].apply(sum)

    # Percentage of ligands that fullfill the druglikeness criteria
    drug_likeness_ratio = round(drug_likeness / recombined_ligands_properties.shape[0] * 100, 0)

    return drug_likeness_ratio
