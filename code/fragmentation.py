from rdkit import Chem
from rdkit.Chem import BRICS

from classes import Fragment
from functions import getDummyLabel


# returns atom numbers of BRICS fragments
def FindBRICSFragments(mol):

    atomTuples = [bond[0] for bond in list(BRICS.FindBRICSBonds(mol))]
    bonds = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atomTuples]
    brokenMol = Chem.FragmentOnBonds(mol, bonds)

    # brokenMol = BRICS.BreakBRICSBonds(mol)
    # Draw.MolToFile(brokenMol, 'test/3w2s_ligand_broken.png')

    fragments = Chem.GetMolFrags(brokenMol)

    return [[atom for atom in fragment if atom < mol.GetNumAtoms()] for fragment in fragments], atomTuples


# given list of atom tuples and a ligand, returns a list of Fragments
# ligand is fragmented at the bonds corresponding to the atom tuples
def getFragmentsFromAtomTuples(atomTuples, BRICSFragments, ligand):

    # get rdkit bonds (NOT BRICS bonds but custom bonds already)
    bonds = [ligand.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atomTuples]
    # fragment ligand at bonds
    fragmentedLigand = Chem.FragmentOnBonds(ligand, bonds)
    # get smiles of fragments
    fragmentSmiles = Chem.MolToSmiles(fragmentedLigand).split('.')
    # get rdkit molecules of fragments
    fragmentMols = Chem.GetMolFrags(fragmentedLigand, asMols=True)  # [Chem.MolFromSmiles(x) for x in fragmentSmiles]
    # get atom numbers (w.r.t. ligand) of fragments
    fragmentAtoms = Chem.GetMolFrags(fragmentedLigand)

    fragments = []

    # iterate over new fragments
    for i, atomNumbers in enumerate(fragmentAtoms):

        # get subpocket corresponding to fragment (Is there a better way?)
        subpocket = [BRICSFragment.subpocket for BRICSFragment in BRICSFragments if atomNumbers[0] in BRICSFragment.atomNumbers][0]
        mol = fragmentMols[i]
        smiles = fragmentSmiles[i]
        # create Fragment object
        fragment = Fragment(mol=mol, smiles=smiles, atomNumbers=atomNumbers, subpocket=subpocket)

        # set atom properties for the created fragment
        for atom in fragment.mol.GetAtoms():
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == '*':
                    atom.SetProp('priority', '2')
                    # getDummyLabel should return the actual atom number w.r.t. the ligand of the dummy atom
                    neighboringSubpocket = [BRICSFragment.subpocket for BRICSFragment in BRICSFragments
                                            if getDummyLabel(neighbor) in BRICSFragment.atomNumbers][0]
                    atom.SetProp('neighboringSubpocket', neighboringSubpocket)

        fragments.append(fragment)

    return fragments
