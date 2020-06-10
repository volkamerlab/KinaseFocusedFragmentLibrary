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


def prepare_chembl(in_file, out_file):

    start = time.time()

    in_file = Path(in_file)
    out_file = Path(out_file)

    logging.basicConfig(filename=out_file.parent / f'{out_file.stem}.log', level=logging.INFO)

    logger.info(f'Read {in_file}...')

    molecules = pd.read_csv(in_file, sep='\t')
    logger.info(f'Number of initial ChEMBL molecules: {molecules.shape[0]}')
    # downloaded ChEMBL file contains data on: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key
    # drop columns not needed
    molecules.drop(['canonical_smiles', 'standard_inchi_key'], axis='columns', inplace=True)

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

    runtime = time.time() - start
    logger.info(f'Time: {runtime}')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--chembl_downloaded_file', type=str, help='file with downloaded ChEMBL data', required=True)
    parser.add_argument('-o', '--chembl_standardized_file', type=str, help='output file for standardized ChEMBL InChIs', required=True)
    args = parser.parse_args()

    # get paths
    chembl_downloaded_file = Path(args.chembl_downloaded_file)
    chembl_standardized_file = Path(args.chembl_standardized_file)

    # standardize chembl
    prepare_chembl(chembl_downloaded_file, chembl_standardized_file)


if __name__ == "__main__":
    main()
