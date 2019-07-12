from rdkit import Chem
import dask.dataframe as dd
import pandas as pd
from rdkit.Chem import PandasTools


def read_inchis(file):

    print('Read', file)

    inchis = pd.read_csv(file, header=None, names=['inchi'], squeeze=True)

    # remove all protonations
    inchis = inchis.str.replace(r'/p\+1', '')
    # print(inchis[inchis.str.find('/p+1') != -1].values)

    print('Number of ChEMBL molecules:', inchis.shape[0])

    return inchis


def read_scaffolds(files):

    scaffolds = pd.DataFrame()

    for file in files:

        print('Read', file)

        scaffolds = scaffolds.append(pd.read_csv(file, sep='\t'))

    print('Number of kinase inhibitor scaffolds:', scaffolds.shape[0])

    PandasTools.AddMoleculeColumnToFrame(scaffolds, 'BM_SMILES', 'mol')

    return scaffolds


def read_original_ligands(frag_dict):

    print('Read original ligands.')

    kinases_pdbs = set()

    for subpocket in frag_dict:

        for frag in frag_dict[subpocket]:
            kinases_pdbs.add((frag.GetProp('kinase'), frag.GetProp('_Name')))

    smiles = []
    inchis = []
    for kinase, pdb in kinases_pdbs:
        f = '/home/paula/Masterarbeit/data/KLIFS_download/HUMAN/' + kinase + '/' + pdb + '/ligand.mol2'
        ligand = Chem.MolFromMol2File(f)
        # ligand.SetProp('complex_pdb', pdb)
        smiles.append(Chem.MolToSmiles(ligand))
        inchis.append(Chem.MolToInchi(ligand))

    print('Number of original ligands :', len(smiles))

    ligands = pd.DataFrame(data=inchis, dtype=str, columns=['inchi'])
    # remove all protonations
    ligands['inchi'] = ligands.inchi.str.replace(r'/p\+1', '')

    # add molecule column
    ligands['smiles'] = smiles
    PandasTools.AddMoleculeColumnToFrame(ligands, 'smiles', 'mol')

    return ligands


def inchi_to_smiles(inchi):

    try:
        smiles = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
    except:
        smiles = None

    return smiles


def read_chembl(file):

    print('Read ChEMBL.')

    out_file = '../../data/chembl/chembl.txt'

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
