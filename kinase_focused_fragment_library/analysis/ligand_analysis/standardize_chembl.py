import pandas as pd
import time
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


def standardize_smiles(input_smiles):

    try:
        smiles = rdMolStandardize.StandardizeSmiles(input_smiles)
    except Exception as e:
        print(e, input_smiles)
        return None

    return smiles


def read_chembl(in_file, out_file):

    print('Read', in_file)

    # chembl_id, canonical_smiles, standard_inchi, standard_inchi_key

    mols = pd.read_csv(in_file, sep='\t').head(100)
    # mols = mols.drop(['canonical_smiles', 'chembl_id', 'standard_inchi_key'], axis='columns')

    # mols['smiles'] = mols.standard_inchi.apply(inchi_to_smiles, meta=('smiles', str))
    print('Number of ChEMBL molecules:', mols.shape[0])

    # TODO: standardize molecules instead of removing p+1 from InChI
    # write to file
    # mols['standard_inchi'].compute().to_csv(out_file, header=0, index=0)

    start = time.time()
    chembl = mols['canonical_smiles'].apply(standardize_smiles)
    print(time.time()-start)

    chembl = chembl.dropna(how='any')

    print('Number of filtered ChEMBL molecules:', len(chembl))

    chembl.to_csv(out_file, header=0, index=0)

    return chembl


read_chembl('chembl_25_chemreps.txt', 'chembl_smiles.txt')