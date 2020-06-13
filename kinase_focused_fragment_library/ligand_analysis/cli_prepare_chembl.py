"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This is the CLI to prepare the ChEMBL dataset (only run once with new ChEMBL dataset from ChEMBL website).
Preparation includes the standardization of InChIs and the removal of ligands that have any missing data.
The standardized dataset is written to a CSV file (ChEMBL compound IDs and the standardized InChIs per ligand).
"""

import argparse
import logging
from pathlib import Path
import time

import pandas as pd
from rdkit.Chem import PandasTools

from .utils import standardize_mol, convert_mol_to_inchi

logger = logging.getLogger(__name__)


def main():
    """
    Main CLI function to standardize a ChEMBL dataset downloaded from ChEMBL and write output to csv file.
    """

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

    # standardize chembl
    prepare_chembl(chembl_downloaded_file, chembl_standardized_file)

    # get script runtime
    runtime = time.time() - start
    logger.info(f'Time: {runtime}')


if __name__ == "__main__":
    main()
