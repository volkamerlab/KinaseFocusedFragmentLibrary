import argparse
import logging
from pathlib import Path
import time

from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)

from kinase_focused_fragment_library.ligand_analysis.utils import read_fragment_library, read_original_ligands, read_chembl_ligands
from kinase_focused_fragment_library.ligand_analysis.base import CombinatorialLibraryAnalyzer

logger = logging.getLogger(__name__)


def main():

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
    parser.add_argument('-chembl', type=str, help='file with standardized chembl data (InChIs)', required=True)
    parser.add_argument('-klifs', type=str, help='path to KLIFS_download folder (original ligands)', required=True)
    parser.add_argument('-o', '--combinatoriallibrary', type=str, help='output path', required=True)
    args = parser.parse_args()

    # set paths
    path_fragment_library = Path(args.fragmentlibrary)
    path_chembl_data = Path(args.chembl)
    path_klifs_data = Path(args.klifs) / 'KLIFS_download'
    path_combinatorial_library = Path(args.combinatoriallibrary)

    # configure logging file
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(path_combinatorial_library / f'combinatorial_library.log',),
            logging.StreamHandler()
        ]
    )

    # get start time of script
    start = time.time()

    # load fragment library
    subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']
    fragment_library = read_fragment_library(path_fragment_library, subpockets)

    # load standardized ChEMBL InChIs
    chembl = read_chembl_ligands(path_chembl_data)

    # load original ligands from KLIFS
    original_ligands = read_original_ligands(fragment_library, path_klifs_data)

    # construct and analyze ligands (write to json file)
    analyzer = CombinatorialLibraryAnalyzer()
    analyzer.run(
        fragment_library,
        original_ligands,
        chembl,
        path_combinatorial_library
    )

    runtime = time.time() - start
    logger.info(f'Time: {runtime}')


if __name__ == "__main__":
    main()
