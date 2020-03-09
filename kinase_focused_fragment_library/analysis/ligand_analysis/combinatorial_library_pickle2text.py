"""
Reformat and save combinatorial library from pickle file to text files.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.append('../../recombination')  # Not pretty but pickle file can only be loaded with explicit paths
from pickle_loader import pickle_loader


def main():
    """
    Reformat and save combinatorial library from pickle file to text files.
    """

    parser = argparse.ArgumentParser(
        description='Reformat and save combinatorial library from pickle file to text files..')
    parser.add_argument('-f', type=str, help='path to combinatorial library pickle file.', required=True)

    args = parser.parse_args()

    properties, fragment_ids, fragment_bonds = get_combinatorial_library_data(
        Path(args.f),
        None
    )

    save_combinatorial_library_data(
        properties,
        fragment_ids,
        fragment_bonds,
        Path(args.f).parent  # save files next to pickle file
    )


def get_combinatorial_library_data(combinatorial_library_file, ligand_ids=None):
    """
    Get ligand properties, ligand fragment IDs and ligand fragment bonds for combinatorial library.

    Parameters
    ----------
    combinatorial_library_file : pathlib.Path
        Path to pickle file containing combinatorial library.
    ligand_ids : None or list of int
        If None, full library is processed (default). Alternatively, pass list of ligand IDs.

    Returns
    -------
    tuple
        - numpy.ndarray: Ligand properties for ligands in combinatorial library.
        - dict of list of str: Ligand fragment IDs (values) for ligands in combinatorial library (keys).
        - dict of list of list of str: Ligand fragment atom pairs describing bonds (value) for ligands in combinatorial 
          library (keys).
    """

    if ligand_ids is None:
        return _get_combinatorial_library_data_all(combinatorial_library_file)
    else:
        return _get_combinatorial_library_data_by_index(combinatorial_library_file, ligand_ids)


def save_combinatorial_library_data(properties, fragment_ids, fragment_bonds, output_path):
    """
    Save ligand properties, ligand fragments IDs and bonds to disc in text files.
    
    Parameters
    ----------
    properties : numpy.ndarray
        Ligand properties for ligands in combinatorial library.
    fragment_ids : dict of list of str
        Ligand fragment IDs (value) for ligands in combinatorial library (keys).
    fragment_bonds : dict of list of list of str
        Ligand fragment atom pairs describing bonds (value) for ligands in combinatorial library (keys).
    output_path : pathlib.Path
        Path to output folder.
    """
    _save_properties(properties, output_path)
    _save_fragments(fragment_ids, fragment_bonds, output_path)


def _get_combinatorial_library_data_all(combinatorial_library_file):
    """
    Get ligand properties, ligand fragment IDs and ligand fragment bonds for combinatorial library.

    Parameters
    ----------
    combinatorial_library_file : pathlib.Path
        Path to pickle file containing combinatorial library.

    Returns
    -------
    tuple
        - numpy.ndarray: Ligand properties for ligands in combinatorial library.
        - dict of list of str: Ligand fragment IDs (values) for ligands in combinatorial library (keys).
        - dict of list of list of str: Ligand fragment atom pairs describing bonds (value) for ligands in combinatorial 
          library (keys).
    """

    combinatorial_library_file = Path(combinatorial_library_file)
    
    properties = np.empty((0, 11), int)
    fragment_ids = {}
    fragment_bonds = {}

    with open(combinatorial_library_file, 'rb') as pickle_file:

        for i, ligand in enumerate(pickle_loader(pickle_file)):

            if i % 1000000 == 0:
                print(i)

            # get ligand properties
            ligand_properties = _get_ligand_properties(i + 1, ligand)

            properties = np.append(
                properties,
                np.array([ligand_properties]),
                axis=0
            )

            # get fragment bonds and ids
            fragment_ids[i + 1] = _get_ligand_fragment_ids(ligand)
            fragment_bonds[i + 1] = _get_ligand_fragment_bonds(ligand)

    return properties, fragment_ids, fragment_bonds


def _get_combinatorial_library_data_by_index(combinatorial_library_file, ligand_ids):
    """
    Get ligand properties, ligand fragment IDs and ligand fragment bonds for combinatorial library.

    Parameters
    ----------
    combinatorial_library_file : pathlib.Path
        Path to pickle file containing combinatorial library.
    ligand_ids : None or list of int
        If None, full library is processed (default). Alternatively, pass list of ligand IDs.

    Returns
    -------
    tuple
        - numpy.ndarray: Ligand properties for ligands in combinatorial library.
        - dict of list of str: Ligand fragment IDs (values) for ligands in combinatorial library (keys).
        - dict of list of list of str: Ligand fragment atom pairs describing bonds (value) for ligands in combinatorial 
          library (keys).
    """

    combinatorial_library_file = Path(combinatorial_library_file)
    
    properties = np.empty((0, 11), int)
    fragment_ids = {}
    fragment_bonds = {}

    # iterate over ligand ids of interest
    ligand_ids_iterator = iter(sorted(ligand_ids))
    ligand_id = next(ligand_ids_iterator)

    with open(combinatorial_library_file, 'rb') as pickle_file:

        # iterate over ligands (up to maximal ligand id of interest)
        for i in range(max(ligand_ids) + 1):

            ligand = next(pickle_loader(pickle_file))

            # store ligand corresponding to entry in list of ligands of interest
            if i + 1 == ligand_id:

                # get ligand properties
                ligand_properties = _get_ligand_properties(i + 1, ligand)

                properties = np.append(
                    properties,
                    np.array([ligand_properties]),
                    axis=0
                )

                # get fragment bonds and ids
                fragment_ids[i + 1] = _get_ligand_fragment_ids(ligand)
                fragment_bonds[i + 1] = _get_ligand_fragment_bonds(ligand)

                # get next ligand index of interest
                if ligand_id < max(ligand_ids):
                    ligand_id = next(ligand_ids_iterator)

    return properties, fragment_ids, fragment_bonds


def _get_ligand_properties(ligand_id, ligand):
    """
    Get ligand properties from ligand object.
    
    Parameters
    ----------
    ligand_id : int
        Ligand ID.
    ligand : kinase_focused_fragment_library.analysis.ligand_analysis.Result.Result
        Ligand object.

    Returns
    -------
    numpy.ndarray
        Ligand properties for ligands in combinatorial library.
    """

    return [
        ligand_id,
        ligand.mwt,
        ligand.hba,
        ligand.hbd,
        ligand.logp,
        ligand.lipinski,
        ligand.n_atoms,
        ligand.n_subpockets,
        ligand.original,
        ligand.original_sub,
        ligand.chembl_match
    ]


def _get_ligand_fragment_ids(ligand):
    """
    Get ligand fragment IDs from ligand object.

    Parameters
    ----------
    ligand : kinase_focused_fragment_library.analysis.ligand_analysis.Result.Result
        Ligand object.

    Returns
    -------
    dict of list of str
        Ligand fragment IDs (values) for ligands in combinatorial library (keys).
    """

    return list(ligand.meta.frag_ids)


def _get_ligand_fragment_bonds(ligand):
    """
    Get ligand fragment bonds from ligand object.

    Parameters
    ----------
    ligand : kinase_focused_fragment_library.analysis.ligand_analysis.Result.Result
        Ligand object.

    Returns
    -------
    dict of list of list of str
    Ligand fragment atom pairs describing bonds (value) for ligands in combinatorial library (keys).
    """

    return [list(i) for i in ligand.meta.bonds]


def _save_properties(properties, output_path):
    """
    Save ligand properties to disc in text files.

    Parameters
    ----------
    properties : numpy.ndarray
        Ligand properties for ligands in combinatorial library.
    output_path : pathlib.Path
        Path to output folder.
    """

    np.savetxt(
        output_path / 'combinatorial_library_properties.csv',
        properties,
        delimiter=',',
        fmt='%i'
    )


def _save_fragments(fragment_ids, fragment_bonds, output_path):
    """
    Save ligand ligand fragments IDs and bonds to disc in text files.

    Parameters
    ----------
    fragment_ids : dict of list of str
        Ligand fragment IDs (value) for ligands in combinatorial library (keys).
    fragment_bonds : dict of list of list of str
        Ligand fragment atom pairs describing bonds (value) for ligands in combinatorial library (keys).
    output_path : pathlib.Path
        Path to output folder.
    """

    output_path = Path(output_path)

    with open(output_path / 'combinatorial_library_fragment_ids.json', 'w') as f:
        json.dump(fragment_ids, f)

    with open(output_path / 'combinatorial_library_fragment_bonds.json', 'w') as f:
        json.dump(fragment_bonds, f)


if __name__ == "__main__":
    main()
