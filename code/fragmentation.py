from rdkit import Chem
from rdkit.Chem import BRICS

from classes import Fragment


# returns atom numbers of BRICS fragments + bond tuples
def FindBRICSFragments(mol):

    atomTuples = [bond[0] for bond in list(BRICS.FindBRICSBonds(mol))]
    bonds = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atomTuples]
    brokenMol = Chem.FragmentOnBonds(mol, bonds, addDummies=False)

    # brokenMol = BRICS.BreakBRICSBonds(mol)
    # Draw.MolToFile(brokenMol, 'test/3w2s_ligand_broken.png')

    fragments = Chem.GetMolFrags(brokenMol)

    return fragments, atomTuples


# given a list of atom tuples and a ligand, returns a list of Fragment objects
# ligand is fragmented at the bonds corresponding to the atom tuples
def getFragmentsFromAtomTuples(atomTuples, BRICSFragments, ligand):

    # get rdkit bonds (NOT BRICS bonds but custom bonds already)
    bonds = [ligand.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atomTuples]
    # fragment ligand at bonds and keep dummy atoms
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
        for a, atom in enumerate(fragment.mol.GetAtoms()):

            # if atom is not a dummy atom
            if atom.GetSymbol() != '*':
                # set atom number within the entire molecule as property of the fragment atom
                # IS THIS ALWAYS TRUE? (Does order of atoms always stay the same after fragmentation?)
                atom.SetProp('atomNumber', str(fragment.atomNumbers[a]))
                atom.SetProp('neighboringSubpocket', 'None')
                atom.SetProp('priority', '1')

            # if atom = dummy atom
            else:
                atom.SetProp('priority', '0')
                atom.SetProp('neighboringSubpocket', 'None')
                atom.SetProp('subpocket', 'None')
                # neighbor = atom next to a bond (Can several neighbors exist?)
                for neighbor in atom.GetNeighbors():
                    neighbor.SetProp('priority', '2')

                neighborAtom = int(neighbor.GetProp('atomNumber'))
                # get and set atom number w.r.t ligand of the dummy atom
                bondAtoms = [atomTuple for atomTuple in atomTuples if neighborAtom in atomTuple][0]
                dummyAtom = [atomNumber for atomNumber in bondAtoms if atomNumber != neighborAtom][0]
                atom.SetProp('atomNumber', str(dummyAtom))
                # get and set neighboring subpocket of the dummy atom
                neighboringSubpocket = [BRICSFragment.subpocket for BRICSFragment in BRICSFragments
                                        if dummyAtom in BRICSFragment.atomNumbers][0]
                for neighbor in atom.GetNeighbors():
                    neighbor.SetProp('neighboringSubpocket', neighboringSubpocket)

        fragments.append(fragment)

    return fragments
