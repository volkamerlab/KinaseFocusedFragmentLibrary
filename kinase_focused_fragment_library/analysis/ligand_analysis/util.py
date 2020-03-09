"""
Utility functions to work with the combinatorial library.
"""

from pathlib import Path

import numpy as np
import pandas as pd


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
