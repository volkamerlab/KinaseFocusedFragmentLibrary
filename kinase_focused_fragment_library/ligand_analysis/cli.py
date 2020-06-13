"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

Command line interfaces to (i) prepare the ChEMBL dataset and (ii) analyze the combinatorial library.
"""

import argparse
import logging
from pathlib import Path
import time

from .base import ChemblPreparer, CombinatorialLibraryAnalyzer
from .utils import read_fragment_library, read_original_ligands, read_chembl_ligands

logger = logging.getLogger(__name__)


def prepare_chembl():
    """
    Command line interface to prepare ChEMBL dataset for the analysis of the combinatorial library.
    """

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--chembl_downloaded_file', type=str, help='file with downloaded ChEMBL data', required=True)
    parser.add_argument('-o', '--chembl_standardized_file', type=str, help='output file for standardized ChEMBL InChIs', required=True)
    args = parser.parse_args()

    # get paths
    chembl_downloaded_file = Path(args.chembl_downloaded_file)
    chembl_standardized_file = Path(args.chembl_standardized_file)

    # configure logging file
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(chembl_standardized_file.parent / f'{chembl_standardized_file.stem}.log'),
            logging.StreamHandler()
        ]
    )

    # get start time of script
    start = time.time()

    # prepare ChEMBL dataset
    chembl_preparer = ChemblPreparer()
    chembl_preparer.run(chembl_downloaded_file, chembl_standardized_file)

    # get script runtime
    runtime = time.time() - start
    logger.info(f'Time: {runtime}')


def analyze_combinatorial_library():
    """
    Command line interface to analyse the combinatorial library.
    """

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
    parser.add_argument('-chembl', type=str, help='file with standardized chembl data (InChIs)', required=True)
    parser.add_argument('-klifs', type=str, help='path to KLIFS_download folder (original ligands)', required=True)
    parser.add_argument('-o', '--combinatoriallibrary', type=str, help='output path', required=True)
    parser.add_argument('-c', '--ncores', type=str, help='number of cores', required=True)
    args = parser.parse_args()

    # set paths and number of cores
    path_fragment_library = Path(args.fragmentlibrary)
    path_chembl_data = Path(args.chembl)
    path_klifs_data = Path(args.klifs) / 'KLIFS_download'
    path_combinatorial_library = Path(args.combinatoriallibrary)
    n_cores = int(args.ncores)

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
        path_combinatorial_library,
        n_cores
    )

    runtime = time.time() - start
    logger.info(f'Time: {runtime}')
