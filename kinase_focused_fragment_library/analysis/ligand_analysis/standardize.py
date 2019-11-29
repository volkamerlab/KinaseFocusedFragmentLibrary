from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
import argparse

RDLogger.DisableLog('rdApp.*')


# TODO: include try-excepts here!
def standardize_mol(mol):
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
    mol = rdMolStandardize.Normalize(mol)
    mol = rdMolStandardize.Reionize(mol)
    u = rdMolStandardize.Uncharger()
    mol = u.uncharge(mol)
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
