import argparse
import pandas as pd

from .utils import standardize_inchi


def prepare_chembl(in_file, out_file):

    print('Read', in_file)

    # chembl_id, canonical_smiles, standard_inchi, standard_inchi_key

    mols = pd.read_csv(in_file, sep='\t')
    chembl = mols.drop(['canonical_smiles', 'chembl_id', 'standard_inchi_key'], axis='columns')

    print('Number of ChEMBL molecules:', mols.shape[0])

    chembl['standard_inchi_new'] = chembl['standard_inchi'].apply(standardize_inchi)
    chembl['diff'] = chembl.apply(lambda x: x['standard_inchi'] != x['standard_inchi_new'], axis=1)
    print('Standardized ChEMBL molecules:', sum(chembl['diff']))

    chembl = chembl['standard_inchi_new']
    chembl = chembl.dropna(how='any')

    chembl.to_csv(out_file, header=0, index=0)

    print('Number of filtered ChEMBL molecules:', len(chembl), mols.shape[0]-len(chembl))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--chembl_downloaded_file', type=str, help='file with downloaded ChEMBL data', required=True)
    parser.add_argument('-o', '--chembl_standardized_file', type=str, help='output file for standardized ChEMBL InChIs', required=True)
    args = parser.parse_args()

    # standardize chembl
    prepare_chembl(args.f, args.o)


if __name__ == "__main__":
    main()
