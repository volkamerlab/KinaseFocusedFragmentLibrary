from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw


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
