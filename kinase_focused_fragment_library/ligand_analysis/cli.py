import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)

from kinase_focused_fragment_library.ligand_analysis.utils import read_fragment_library, read_original_ligands
from kinase_focused_fragment_library.ligand_analysis.analyze_results import get_ligands_analysis


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

    # write ligand analysis to json file
    get_ligands_analysis(fragment_library, original_ligands, chembl, in_paths, combinatorial_library_file)


if __name__ == "__main__":
    main()
