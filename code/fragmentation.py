from rdkit import Chem
from rdkit.Chem import BRICS

from classes import Fragment


# returns atom numbers of BRICS fragments + bond tuples
def find_brics_fragments(mol):

    atom_tuples = [bond[0] for bond in BRICS.FindBRICSBonds(mol)]
    # if mol was not fragmented:
    if len(atom_tuples) == 0:
        fragments = [Fragment(atomNumbers=range(mol.GetNumAtoms()), mol=mol)]
        return fragments, atom_tuples
    # else:
    bonds = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atom_tuples]
    broken_mol = Chem.FragmentOnBonds(mol, bonds, addDummies=False)

    fragment_atoms = Chem.GetMolFrags(broken_mol)
    fragment_mols = Chem.GetMolFrags(broken_mol, asMols=True)

    fragments = [Fragment(atomNumbers=n, mol=m) for (n, m) in zip(fragment_atoms, fragment_mols)]

    return fragments, atom_tuples


# given a list of atom tuples, BRICS fragments, and the ligand, returns a list of Fragment objects
# ligand is fragmented at the bonds corresponding to the atom tuples
def fragmentation(ligand, atom_tuples, brics_fragments):

    # get rdkit bonds (NOT BRICS bonds but custom bonds already)
    bonds = [ligand.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atom_tuples]
    if len(bonds) > 0:
        # fragment ligand at bonds and keep dummy atoms
        fragmented_ligand = Chem.FragmentOnBonds(ligand, bonds)
    else:
        fragmented_ligand = ligand
    # get smiles of fragments
    fragment_smiles = Chem.MolToSmiles(fragmented_ligand).split('.')
    # get rdkit molecules of fragments
    fragment_mols = Chem.GetMolFrags(fragmented_ligand, asMols=True)  # [Chem.MolFromSmiles(x) for x in fragment_smiles]
    # get atom numbers (w.r.t. ligand) of fragments
    fragment_atoms = Chem.GetMolFrags(fragmented_ligand)

    fragments = []

    # iterate over new fragments
    for i, atomNumbers in enumerate(fragment_atoms):

        # get subpocket corresponding to fragment (Is there a better way?)
        subpocket = [brics_fragment.subpocket for brics_fragment in brics_fragments if atomNumbers[0] in brics_fragment.atomNumbers][0]
        # create Fragment object
        fragment = Fragment(mol=fragment_mols[i], smiles=fragment_smiles[i], atomNumbers=atomNumbers, subpocket=subpocket)

        # set atom properties for the created fragment
        for a, atom in enumerate(fragment.mol.GetAtoms()):

            # if atom is not a dummy atom
            if atom.GetSymbol() != '*':
                # set atom number within the entire molecule as property of the fragment atom
                # IS THIS ALWAYS TRUE? (Does order of atoms always stay the same after fragmentation?)
                atom.SetProp('atomNumber', str(fragment.atomNumbers[a]))
                #atom.SetProp('neighboringSubpocket', 'None')
                atom.SetProp('priority', '1')

            # if atom = dummy atom
            else:
                atom.SetProp('priority', '0')
                #atom.SetProp('neighboringSubpocket', 'None')
                # atom.SetProp('subpocket', 'None')
                # neighbor = atom next to a bond (Can several neighbors exist?)
                for neighbor in atom.GetNeighbors():
                    neighbor.SetProp('priority', '2')

                # Problem: Why does this work? Is atomNumber prop defined somewhere else except in above if statement?
                neighbor_atom = int(neighbor.GetProp('atomNumber'))
                # get and set atom number w.r.t ligand of the dummy atom
                bond_atoms = [atomTuple for atomTuple in atom_tuples if neighbor_atom in atomTuple][0]
                dummy_atom = [atomNumber for atomNumber in bond_atoms if atomNumber != neighbor_atom][0]
                atom.SetProp('atomNumber', str(dummy_atom))
                # get and set neighboring subpocket of the dummy atom
                neighboring_subpocket = [BRICSFragment.subpocket for BRICSFragment in brics_fragments
                                         if dummy_atom in BRICSFragment.atomNumbers][0]
                for neighbor in atom.GetNeighbors():
                    neighbor.SetProp('neighboringSubpocket', neighboring_subpocket)

        fragments.append(fragment)

    return fragments
