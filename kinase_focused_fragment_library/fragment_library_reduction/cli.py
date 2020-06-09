import argparse
from pathlib import Path

from rdkit import Chem

from kinase_focused_fragment_library.utils import read_fragment_library


def main():

    # ============================= COMMAND LINE ARGUMENTS ===================================

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
    parser.add_argument('-o', '--fragmentlibraryreduced', type=str, help='output path for reduced fragment library', required=True)
    args = parser.parse_args()

    # ============================= PATHS =====================================================

    PATH_FRAG_LIB = Path(args.fragmentlibrary)
    PATH_FRAG_LIB_REDUCED = Path(args.fragmentlibraryreduced)
    PATH_FRAG_LIB_REDUCED.mkdir(parents=True, exist_ok=True)

    # ============================= INPUT DATA ================================================

    fragment_library = read_fragment_library(PATH_FRAG_LIB)

    # ============================= OUTPUT DATA kff===============================================

    _save_top_n_fragments(fragment_library, 5, PATH_FRAG_LIB_REDUCED)


def _save_top_n_fragments(fragment_library, n, path_fragment_library_reduced):

    for subpocket, fragments in fragment_library.items():

        if subpocket != 'X':

            with open(path_fragment_library_reduced / f'{subpocket}.sdf', 'w') as f:

                w = Chem.SDWriter(f)
                for fragment in fragments.ROMol_original[:n]:
                    w.write(fragment)
                w.close()
