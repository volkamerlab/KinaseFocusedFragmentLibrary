from rdkit import Chem
import dask.dataframe as dd
import pandas as pd


def read_mols(file):

    print('Read', file)

    mols = pd.read_csv(file, header=None, names=['string'])

    # print(smiles[smiles.smiles.str.find('+') != -1])

    print('Number of ChEMBL molecules:', mols.shape[0])

    return mols


def inchi_to_smiles(inchi):

    try:
        smiles = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
    except:
        smiles = None

    return smiles


def read_chembl(file):

    print('Read ChEMBL.')

    out_file = '../../chembl/chembl.txt'

    # chembl_id, canonical_smiles, standard_inchi, standard_inchi_key

    mols = dd.read_csv(file, sep='\t')
    # convert to canonical rdkit smiles
    mols = mols.drop(['canonical_smiles', 'chembl_id', 'standard_inchi_key'], axis='columns')

    # mols['smiles'] = mols.standard_inchi.apply(inchi_to_smiles, meta=('smiles', str))
    print('Number of ChEMBL molecules:', mols['standard_inchi'].compute().shape[0])
    mols = mols.dropna(how='any')
    print('Number of filtered ChEMBL molecules:', mols['standard_inchi'].compute().shape[0])

    # write to file
    mols['standard_inchi'].compute().to_csv(out_file, header=0, index=0)

    return out_file


# print(read_chembl('/home/paula/Downloads/chembl_25_chemreps.txt'))
