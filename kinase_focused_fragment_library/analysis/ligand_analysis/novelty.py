from rdkit import Chem
import pandas as pd

from .standardize import standardize_mol, standardize_inchi


def read_chembl(in_file):

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

    chembl.to_csv('chembl_standardized_inchi', header=0, index=0)

    print('Number of filtered ChEMBL molecules:', len(chembl), mols.shape[0]-len(chembl))

    return chembl


def read_original_ligands(frag_dict, path_to_klifs):

    print('Read original ligands.')

    kinases_pdbs = set()

    for subpocket in frag_dict:

        for frag in frag_dict[subpocket]:
            kinases_pdbs.add((frag.GetProp('kinase'), frag.GetProp('_Name')))

    inchis = []
    mols = []
    for kinase, pdb in kinases_pdbs:
        f = path_to_klifs / ('HUMAN/' + kinase + '/' + pdb + '/ligand.mol2')
        ligand = Chem.MolFromMol2File(str(f))

        # standardization
        ligand = standardize_mol(ligand)
        # if ligand could not be standardized, skip
        if not ligand:
            print('Ligand could not be standardized: ', pdb)
            return

        mols.append(ligand)
        inchi = Chem.MolToInchi(ligand)
        inchis.append(inchi)

    print('Number of original ligands :', len(inchis))

    ligands = pd.DataFrame(data=inchis, dtype=str, columns=['inchi'])
    # add molecule column
    ligands['mol'] = mols

    return ligands
