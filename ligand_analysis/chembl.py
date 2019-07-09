from rdkit import Chem
import dask.dataframe as dd


def read_chembl(file):

    print('Read ChEMBL.')

    # chembl_id, canonical_smiles, standard_inchi, standard_inchi_key

    mols = dd.read_csv(file, sep='\t')
    # convert to canonical rdkit smiles
    mols = mols.drop(['canonical_smiles', 'chembl_id', 'standard_inchi_key'], axis='columns')

    print('Number of ChEMBL molecules: ', mols.compute().shape[0])

    return mols


# import apsw
#
#
# # the extension is usually loaded right after the connection to the
# # database
# connection = apsw.Connection('chembldb.sql')
# connection.enableloadextension(True)
# connection.loadextension('~/chemicalite/')
# connection.enableloadextension(False)
