import argparse
import multiprocessing as mp
from pathlib import Path
import time
import pickle

# import matplotlib.pyplot as plt
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
    fragments = read_fragment_library(Path(args.fragmentlibrary), subpockets)

    # read standardized chembl inchis
    print('Read', args.chembl)
    chembl = pd.read_csv(args.chembl, header=None, names=['standard_inchi'])
    print('Number of ChEMBL molecules:', chembl.shape[0])

    # original ligands from KLIFS
    path_to_klifs = Path(args.klifs) / 'KLIFS_download'
    original_ligands = read_original_ligands(fragments, path_to_klifs)

    # output file
    combinatorial_library_folder = Path(args.combinatoriallibrary)
    combinatorial_library_file = combinatorial_library_folder / 'combinatorial_library.pickle'

    # objects create by the recombination algorithm
    path_to_results = combinatorial_library_folder / 'results'
    in_paths = list(path_to_results.glob('*.pickle'))

    # ================================ INITIALIZE =========================================


    n_per_sp, n_filtered_per_sp = {}, {}
    for subpocket in subpockets:
        n_per_sp[subpocket] = 0
        n_filtered_per_sp[subpocket] = 0

    n_sp, n_filtered_sp = {}, {}
    for i in range(len(subpockets)):
        n_sp[i+1] = 0
        n_filtered_sp[i+1] = 0

    results = []

    n_processes = 2  # mp.cpu_count()
    print("Number of processors: ", n_processes)
    pool = mp.Pool(n_processes)

    # ========================= CONSTRUCT AND ANALYZE LIGANDS ==============================

    start = time.time()

    combinatorial_library_file = combinatorial_library_file.open('wb')

    # iterate over ligands
    for in_path in in_paths:

        print(str(in_path))
        with open(str(in_path), 'rb') as pickle_in:

            tmp_results = pool.starmap(analyze_result, [(meta, fragments, original_ligands, chembl) for meta in pickle_loader(pickle_in)])

            # store in combinatorial library
            for result in tmp_results:
                pickle.dump(result, combinatorial_library_file)

            combinatorial_library_file.flush()

            results.extend(tmp_results)

    combinatorial_library_file.close()

    runtime = time.time() - start
    print('Time: ', runtime)


if __name__ == "__main__":
    main()
