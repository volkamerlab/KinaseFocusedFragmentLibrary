from rdkit import Chem
from rdkit.Chem import BRICS

from classes import Fragment


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
def GetFragmentsFromAtomTuples(atomTuples, BRICSFragments, ligand):

    # get rdkit bonds (NOT BRICS bonds but custom bonds already)
    bonds = [ligand.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atomTuples]
    # get subpocket labels for dummy atoms
    labels = []
    for i, bond in enumerate(bonds):
        beginAtom = atomTuples[i][0]
        endAtom = atomTuples[i][1]
        # set neighboring subpockets as labels for dummy atoms
        firstSubpocket = [fragment.subpocket for fragment in BRICSFragments if endAtom in fragment.atomNumbers][0]
        secondSubpocket = [fragment.subpocket for fragment in BRICSFragments if beginAtom in fragment.atomNumbers][0]
        labels.append((int(firstSubpocket), int(secondSubpocket)))

    # fragment ligand at bonds
    fragmentedLigand = Chem.FragmentOnBonds(ligand, bonds, dummyLabels=labels)
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
        fragments.append(Fragment(mol=mol, smiles=smiles, atomNumbers=atomNumbers, subpocket=subpocket))

    return fragments
