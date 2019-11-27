from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools

from standardize import standardize_smiles


def read_smiles(file):

    print('Read', file)

    smiles = pd.read_csv(file, header=None, names=['smiles'])

    return smiles


def read_scaffolds(files):

    scaffolds = pd.DataFrame()

    for file in files:

        print('Read', file)

        scaffolds = scaffolds.append(pd.read_csv(file, sep='\t'))

    print('Number of kinase inhibitor scaffolds:', scaffolds.shape[0])

    PandasTools.AddMoleculeColumnToFrame(scaffolds, 'BM_SMILES', 'mol')

    return scaffolds


def read_original_ligands(frag_dict, path_to_klifs):

    print('Read original ligands.')

    kinases_pdbs = set()

    for subpocket in frag_dict:

        for frag in frag_dict[subpocket]:
            kinases_pdbs.add((frag.GetProp('kinase'), frag.GetProp('_Name')))

    smiles = []
    # inchis = []
    for kinase, pdb in kinases_pdbs:
        f = path_to_klifs / ('HUMAN/' + kinase + '/' + pdb + '/ligand.mol2')
        ligand = Chem.MolFromMol2File(str(f))
        # standardization
        # ligand = standardize_mol(ligand)
        # ligand.SetProp('complex_pdb', pdb)
        s = Chem.MolToSmiles(ligand)
        s = standardize_smiles(s)
        smiles.append(s)
        # inchis.append(Chem.MolToInchi(ligand))

    print('Number of original ligands :', len(smiles))

    # ligands = pd.DataFrame(data=inchis, dtype=str, columns=['inchi'])

    # TODO: standardize molecules instead of removing p+1 from InChI
    # remove all protonations
    # ligands['inchi'] = ligands.inchi.str.replace(r'/p\+1', '')

    ligands = pd.DataFrame(data=smiles, dtype=str, columns=['smiles'])

    # add molecule column
    # ligands['smiles'] = smiles
    PandasTools.AddMoleculeColumnToFrame(ligands, 'smiles', 'mol')

    return ligands
