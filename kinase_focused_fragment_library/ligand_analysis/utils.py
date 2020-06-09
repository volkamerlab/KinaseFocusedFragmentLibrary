import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog('rdApp.*')


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


def standardize_mol(mol):
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.Normalize(mol)
        mol = rdMolStandardize.Reionize(mol)
        u = rdMolStandardize.Uncharger()
        mol = u.uncharge(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        return mol
    except Exception as e:
        print(e)
        return None


def standardize_inchi(input_inchi):
    try:
        mol = Chem.MolFromInchi(input_inchi)
    except:
        print('ERROR in MolFromInchi:', input_inchi)
        return None
    if not mol:
        print('ERROR in MolFromInchi:', input_inchi)
        return None
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.Normalize(mol)
        mol = rdMolStandardize.Reionize(mol)
        u = rdMolStandardize.Uncharger()
        mol = u.uncharge(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except:
        print('ERROR in standardization:', input_inchi)
        return None
    try:
        inchi = Chem.MolToInchi(mol)
    except:
        print('ERROR in MolToInchi:', input_inchi)
        return None
    return inchi