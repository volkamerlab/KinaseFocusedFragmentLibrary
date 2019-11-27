from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
import argparse

RDLogger.DisableLog('rdApp.*')


def standardize_smiles(input_smiles):

    try:
        smiles = rdMolStandardize.StandardizeSmiles(input_smiles)
    except Exception as e:
        print(e, input_smiles)
        return None

    return smiles


def standardize_mol(mol):
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
    mol = rdMolStandardize.Normalize(mol)
    mol = rdMolStandardize.Reionize(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    return mol


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
    except:
        print('ERROR in SanitizeMol:', input_inchi)
        return None
    mol = Chem.RemoveHs(mol)
    mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
    mol = rdMolStandardize.Normalize(mol)
    mol = rdMolStandardize.Reionize(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    return Chem.MolToInchi(mol)


def command_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("inchi", type=str)
    args = parser.parse_args()
    inchi = args.inchi
    return inchi


if __name__ == "__main__":
    input_inchi = command_parser()
    inchi = standardize_mol(input_inchi)
    print(inchi)
