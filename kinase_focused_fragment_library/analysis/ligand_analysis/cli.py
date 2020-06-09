import argparse
import json
import multiprocessing as mp
from pathlib import Path
import time

import pandas as pd
from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)

from kinase_focused_fragment_library.analysis.ligand_analysis.construct_ligand import read_fragment_library
from kinase_focused_fragment_library.recombination.pickle_loader import pickle_loader
from kinase_focused_fragment_library.analysis.ligand_analysis.analyze_results import analyze_result
from kinase_focused_fragment_library.analysis.ligand_analysis.novelty import read_original_ligands


def main():

    # ============================= COMMAND LINE ARGUMENTS ===================================

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
    parser.add_argument('-chembl', type=str, help='file with standardized chembl data (InChIs)', required=True)
    parser.add_argument('-klifs', type=str, help='path to KLIFS_download folder (original ligands)', required=True)
    parser.add_argument('-o', '--combinatoriallibrary', type=str, help='output path', required=True)
    args = parser.parse_args()

    # ============================= INPUT DATA ================================================

    subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']
    fragment_library = read_fragment_library(Path(args.fragmentlibrary), subpockets)

    # read standardized chembl inchis
    print('Read', args.chembl)
    chembl = pd.read_csv(args.chembl, header=None, names=['standard_inchi'])
    print('Number of ChEMBL molecules:', chembl.shape[0])

    # original ligands from KLIFS
    path_to_klifs = Path(args.klifs) / 'KLIFS_download'
    original_ligands = read_original_ligands(fragment_library, path_to_klifs)

    # output file
    combinatorial_library_folder = Path(args.combinatoriallibrary)
    combinatorial_library_file = combinatorial_library_folder / 'combinatorial_library.json'

    # objects create by the recombination algorithm
    path_to_results = combinatorial_library_folder / 'results'
    in_paths = list(path_to_results.glob('*.pickle'))

    # ========================= CONSTRUCT AND ANALYZE LIGANDS ==============================

    # write ligand analysis to pickle file
    _construct_and_analyze_ligands(fragment_library, original_ligands, chembl, in_paths, combinatorial_library_file)


def _construct_and_analyze_ligands(fragment_library, original_ligands, chembl, in_paths, combinatorial_library_file):
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
                analyze_result,
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


if __name__ == "__main__":
    main()
