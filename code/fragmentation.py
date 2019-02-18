from rdkit import Chem
from rdkit.Chem import BRICS


# -------- Find BRICS fragments ---------

# returns atom numbers of BRICS fragments
def FindBRICSFragments(mol, asMols=False):

    brokenMol = BRICS.BreakBRICSBonds(mol)
    # Draw.MolToFile(brokenMol, 'test/3w2s_ligand_broken.png')

    fragments = Chem.GetMolFrags(brokenMol)

    if asMols:

        fragmentsAsMols = Chem.GetMolFrags(brokenMol, asMols=True)

        return [[atom for atom in fragment if atom < mol.GetNumAtoms()] for fragment in fragments], fragmentsAsMols

    else:

        return [[atom for atom in fragment if atom < mol.GetNumAtoms()] for fragment in fragments]


