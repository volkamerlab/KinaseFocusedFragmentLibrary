from rdkit import Chem
from rdkit.Chem import BRICS

from .classes import Fragment


def find_brics_fragments(mol):

    """
    Carries out the BRICS fragmentation algorithm and returns the BRICS fragments and bonds

    Parameters
    ----------
    mol: RDKit Mol object
        molecule to be fragmented

    Returns
    -------
    fragments: list(Fragment)
        list of the resulting BRICS fragments as Fragment objects
    brics_bonds: list(tuple(tuple(int), tuple(str)))
        list of tuples, where each tuple represents a BRICS bond between two atoms in the ligand
        - first tuple: atom indices
        - second tuple: BRICS environment types
    """

    brics_bonds = list(BRICS.FindBRICSBonds(mol))  # convert generator to list for later purposes
    atom_tuples = [bond[0] for bond in brics_bonds]

    # if mol was not fragmented:
    if len(atom_tuples) == 0:
        fragments = [Fragment(atomNumbers=range(mol.GetNumAtoms()), mol=mol, environment='na')]
        return fragments, atom_tuples
    # else:
    bonds = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atom_tuples]
    broken_mol = Chem.FragmentOnBonds(mol, bonds, addDummies=False)

    fragment_atoms = Chem.GetMolFrags(broken_mol)
    fragment_mols = Chem.GetMolFrags(broken_mol, asMols=True)

    fragments = [Fragment(atomNumbers=n, mol=m) for (n, m) in zip(fragment_atoms, fragment_mols)]

    return fragments, brics_bonds


def fragment_between_atoms(mol, atom_tuples):

    """
    Fragmentation of the given molecule between the given atom tuples

    Parameters
    ----------
    mol: RDKit Mol object
        molecule to be fragmented
    atom_tuples: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand

    Returns
    -------
    fragment_mols: list(Mol)
        list of the resulting fragments as RDKit Mol objects
    fragment_atoms: list(list(int))
        list of the resulting fragments as atom indices w.r.t. the original molecule
    """

    # get bonds corresponding to atom tuples
    bonds = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atom_tuples]
    if len(bonds) > 0:
        # fragment ligand at bonds and keep dummy atoms
        fragmented_ligand = Chem.FragmentOnBonds(mol, bonds)
    else:
        fragmented_ligand = mol
    # get rdkit molecules of fragments
    fragment_mols = Chem.GetMolFrags(fragmented_ligand, asMols=True)
    # get atom numbers (w.r.t. ligand) of fragments
    fragment_atoms = Chem.GetMolFrags(fragmented_ligand)

    return fragment_mols, fragment_atoms


def set_atom_properties(fragments, atom_tuples, brics_fragments):

    """

    Assigns properties to the atoms of each fragment (in place):
    - subpocket -> Neighboring subpocket is stored at the dummy atoms
    - atom number w.r.t. original ligand
    - BRICS environment type

    Assumption: When iterating over the fragment atoms, dummy atoms are last.

    Parameters
    ----------
    fragments: list(Fragment)
        list of fragments that the ligand consists of
    atom_tuples: list(tuple(int))
            list of atom index tuples, where each tuple represents a bond between two atoms in the ligand (final bonds, NOT BRICS bonds!)
    brics_fragments: list(Fragment)
        list of BRICS fragments of the ligand as Fragment objects

    Returns
    -------
    fragment: Fragment
        Fragment object of the given fragment, including its subpocket

    """

    #print(f'Fragment bonds: {atom_tuples}')

    fragments_with_multiple_subpocket_bonds_per_atom = []

    for fragment in fragments:

        #print(f'\nFragment subpocket: {fragment.subpocket.name}')
        #print(f'Fragment atom numbers: {fragment.atomNumbers}')

        # Make copy of atom tuples (=bonds between atoms which are in different subpockets)
        # Why? We need this later when we iterate over the dummy atoms; each bond, which could be linked to a dummy
        # atom, will be removed. This ensures that always all bonds are considered, even in case there are multiple
        # dummy atoms (and therefore subpocket-bonds for one atom).
        atom_tuples_tmp = atom_tuples.copy()

        # set atom properties for the created fragment
        for atom, atomNumber in zip(fragment.mol.GetAtoms(), fragment.atomNumbers):

            # if atom is not a dummy atom
            if atom.GetSymbol() != '*':

                # get environment type of the brics fragment that the current atom belongs to
                env_type = next(brics_fragment.environment for brics_fragment in brics_fragments if atomNumber in brics_fragment.atomNumbers)

                # set atom number w.r.t. ligand as property of the fragment atom
                # IS THIS ALWAYS TRUE? (Does order of atoms always stay the same after fragmentation?)
                atom.SetIntProp('atomNumber', atomNumber)
                atom.SetProp('subpocket', fragment.subpocket.name)
                atom.SetProp('environment', env_type)

            # if atom = dummy atom
            # the else loop will be entered at the end (dummy atoms are listed at the end)
            else:

                #print(f'> Dummy atom:')
                #print(f'Dummy atom (initially): {atomNumber}')

                # neighbor = atom next to a dummy
                # can several neighbors exist?
                # DS: Since BRICS does not cut in rings, this should not be possible.
                # DS: Added here test in case there are multiple neighbors nevertheless
                neighbor = atom.GetNeighbors()
                if len(neighbor) != 1:
                    raise ValueError(f'Unexpected number of dummy atom neighbors: {len(neighbor)}')
                else:
                    neighbor = neighbor[0]

                # -> This works only because dummy atoms are always last in the iteration
                # DS: The sentence in the line above refers to which code line?

                neighbor_atom = neighbor.GetIntProp('atomNumber')
                #print(f'Atom connected to dummy atom: {neighbor_atom}')

                # get and set atom number w.r.t ligand of the dummy atom
                # (dummy atoms have just been consecutively numbered as additional atoms,
                # we want to have the atom number of the atom that the dummy atom corresponds to)

                #bond_atoms = next(atomTuple for atomTuple in atom_tuples if neighbor_atom in atomTuple)
                # DS: Why do you use next() here? Why not a simple list comprehension?
                # DS: This is problematic when one atom connects multiple times to neighboring subpockets
                # Alternative:
                bond_atoms = [atomTuple for atomTuple in atom_tuples_tmp if neighbor_atom in atomTuple]
                #print(f'Bond(s) involving this dummy atom: {bond_atoms}')

                # Take always first hit (one atom can have one or multiple dummy atoms)
                bond_atoms_first_hit = bond_atoms[0]
                # Trick is: First hit will always be removed from list, thus if there are multiple hits,
                # the other hits will be found in future iterations
                atom_tuples_tmp.remove(bond_atoms_first_hit)

                #dummy_atom = next(atomNumber for atomNumber in bond_atoms if atomNumber != neighbor_atom)
                dummy_atom = [
                    atomNumber for atomNumber in bond_atoms_first_hit if atomNumber != neighbor_atom
                ]
                if len(dummy_atom) != 1:
                    raise ValueError(f'Unexpected number of dummy atoms: {len(dummy_atom)}')
                else:
                    dummy_atom = dummy_atom[0]
                atom.SetIntProp('atomNumber', dummy_atom)
                #print(f'Dummy atom (finally): {dummy_atom}')

                # get and set subpocket of the dummy atom
                #neighboring_subpocket = next(f.subpocket for f in fragments if dummy_atom in f.atomNumbers)
                neighboring_subpocket = [
                    fragment.subpocket
                    for fragment
                    in fragments
                    if dummy_atom in fragment.atomNumbers
                ]
                if len(neighboring_subpocket) != 1:
                    raise ValueError(f'Unexpected number of subpockets: {len(neighboring_subpocket)}')
                else:
                    neighboring_subpocket = neighboring_subpocket[0]
                atom.SetProp('subpocket', neighboring_subpocket.name)
                #print(f'Neighboring subpocket(s): {neighboring_subpocket.name}')

                # dummy atoms do not get an environment type assigned
                atom.SetProp('environment', 'na')

                # Report fragments with mulitple subpocket-connecting-bonds per atom
                if len(bond_atoms) > 1:
                    fragments_with_multiple_subpocket_bonds_per_atom.append(
                        {
                            'dummy_atom': dummy_atom,
                            'subpocket_dummy': neighboring_subpocket.name,
                            'subpocket_fragment': fragment.subpocket.name,
                            'bonds': bond_atoms
                        }
                    )

    return fragments_with_multiple_subpocket_bonds_per_atom
