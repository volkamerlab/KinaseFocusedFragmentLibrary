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

from .utils import standardize_inchi

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
        filename=chembl_standardized_file.parent / f'{chembl_standardized_file.stem}.log',
        level=logging.INFO
    )

    # get start time of script
    start = time.time()

    # standardize chembl
    prepare_chembl(chembl_downloaded_file, chembl_standardized_file)

    # get script runtime
    runtime = time.time() - start
    logger.info(f'Time: {runtime}')


def prepare_chembl(in_file, out_file):
    """
    Prepare ChEMBL dataset, i.e. load data, standardize InChIs, drop entries with any missing information,
    and write data (ChEMBL compound ID, standardized InChI) to output csv file.
    Write logging output to log file with same file name.

    Parameters
    ----------
    in_file : pathlib.Path
        Path to downloaded ChEMBL file.
    out_file : pathlib.Path
        Path to standardized output ChEMBL dataset file.
    """

    in_file = Path(in_file)
    out_file = Path(out_file)

    # read raw ChEMBL data
    logger.info(f'Read {in_file}...')
    molecules = _read_chembl_raw(in_file)
    logger.info(f'Number of initial ChEMBL molecules: {molecules.shape[0]}')

    # standardize InChIs
    logger.info(f'Standardize InChIs...')
    molecules['standard_inchi_new'] = molecules['standard_inchi'].apply(standardize_inchi)

    # check how many InChIs are changed after standardization
    molecules['diff'] = molecules.apply(lambda x: x['standard_inchi'] != x['standard_inchi_new'], axis=1)
    logger.info(f'Number of ChEMBL molecules with changed InChIs after standardization: {sum(molecules["diff"])}')

    # drop "old" standardized molecules and diff column
    molecules.drop(['standard_inchi', 'diff'], axis='columns', inplace=True)

    # rename column
    molecules.rename(columns={'standard_inchi_new': 'standard_inchi'}, inplace=True)

    # drop rows with any data missing
    molecules.dropna(how='any', inplace=True)
    logger.info(f'Number of filtered ChEMBL molecules: {molecules.shape[0]}')

    # save data to file
    logger.info(f'Save to {out_file}...')
    molecules.to_csv(out_file, index=0)


def _read_chembl_raw(path_to_chembl):
    """
    Read the raw ChEMBL data and drop keep only ChEMBL ID and canonical SMILES columns.

    Parameters
    ----------
    path_to_chembl : pathlib.Path
        Path to raw ChEMBL data.

    Returns
    -------
    pandas.DataFrame
        Raw ChEMBL data with ChEMBL ID (chembl_id) and canonical SMILES (canonical_smiles) columns.
    """

    # read data with columns: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key
    molecules = pd.read_csv(path_to_chembl, sep='\t')

    # drop unneeded columns
    molecules.drop(['standard_inchi', 'standard_inchi_key'], axis='columns', inplace=True)

    return molecules


if __name__ == "__main__":
    main()
