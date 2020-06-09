"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This is the CLI to prepare the ChEMBL dataset (only run once with new ChEMBL dataset from ChEMBL website).
Preparation includes the standardization of InChIs and the removal of ligands that have any missing data.
The standardized dataset is written to a CSV file (ChEMBL compound IDs and the standardized InChIs per ligand).
"""

import argparse
import pandas as pd

from .utils import standardize_inchi


def prepare_chembl(in_file, out_file):

    print('Read', in_file)

    mols = pd.read_csv(in_file, sep='\t')
    print('Number of ChEMBL molecules:', mols.shape[0])
    # downloaded ChEMBL file contains data on: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key
    # drop columns not needed
    mols.drop(['canonical_smiles', 'standard_inchi_key'], axis='columns', inplace=True)

    # standardize InChIs
    mols['standard_inchi_new'] = mols['standard_inchi'].apply(standardize_inchi)

    # check how many InChIs are changed after standardization
    mols['diff'] = mols.apply(lambda x: x['standard_inchi'] != x['standard_inchi_new'], axis=1)
    print('ChEMBL molecules with InChIs differing before and after standardization:', sum(mols['diff']))

    # drop "old" standardized molecules and diff column
    mols.drop(['standard_inchi', 'diff'], axis='columns', inplace=True)

    # rename column
    mols.rename(columns={'standard_inchi_new': 'standard_inchi'}, inplace=True)

    # drop rows with any data missing
    chembl = mols.dropna(how='any')
    print(chembl.columns)

    # save data to file
    chembl.to_csv(out_file, index=0)

    print(f'Number of filtered ChEMBL molecules: {len(chembl)}')
    print(f'Number of dropped ChEMBL molecules: {mols.shape[0]-len(chembl)}')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--chembl_downloaded_file', type=str, help='file with downloaded ChEMBL data', required=True)
    parser.add_argument('-o', '--chembl_standardized_file', type=str, help='output file for standardized ChEMBL InChIs', required=True)
    args = parser.parse_args()

    # standardize chembl
    prepare_chembl(args.chembl_downloaded_file, args.chembl_standardized_file)


if __name__ == "__main__":
    main()
