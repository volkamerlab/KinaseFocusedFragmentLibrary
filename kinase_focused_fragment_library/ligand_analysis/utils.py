from rdkit import Chem
import pandas as pd

from .standardize import standardize_mol


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
