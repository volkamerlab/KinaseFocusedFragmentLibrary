import numpy as np
from rdkit import Chem
from rdkit.Chem import BRICS

from kinase_focused_fragment_library.fragmentation.classes import Fragment


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
    Assign properties to the atoms of each fragment (in place):
    - subpocket (in case of dummy atoms: neighboring subpocket is stored)
    - atom number w.r.t. the original ligand
    - BRICS environment type

    Assumption: When iterating over the fragment atoms, dummy atoms are last.

    Parameters
    ----------
    fragments: list(Fragment)
        list of fragments that the ligand consists of
    atom_tuples: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand
        (final bonds, NOT BRICS bonds!)
    brics_fragments: list(Fragment)
        list of BRICS fragments of the ligand as Fragment objects
    """

    for fragment in fragments:

        atoms = fragment.mol.GetAtoms()
        atom_numbers = fragment.atomNumbers

        for atom, atom_number in zip(atoms, atom_numbers):

            if atom.GetSymbol() != '*':  # regular atom

                atom_number, subpocket, environment = _get_atom_properties_regular(
                    atom_number,
                    fragment,
                    brics_fragments
                )
                
                if environment != 'na':

                    if atom_number in [atom_id for atom_id, _ in environment]:
                        environment = next(env for atom_id, env in environment if atom_id == atom_number)
                    else: # not an atom next to a dummy => set environment randomly
                        environment = next(iter(environment))[1]

                atom.SetIntProp('atomNumber', atom_number)
                atom.SetProp('subpocket', subpocket)
                atom.SetProp('environment', environment)

            else:  # dummy atom

                atom_number, subpocket, environment = _get_atom_properties_dummy(
                    atom_number,
                    atom,
                    atom_tuples,
                    fragments
                )
                atom.SetIntProp('atomNumber', atom_number)
                atom.SetProp('subpocket', subpocket)
                atom.SetProp('environment', environment)


def _get_atom_properties_regular(atom_number, fragment, brics_fragments):

    """
    Set atom properties for an atom which is no dummy atom (regular atom):
    each atom's number, subpocket and BRICS environment.

    Parameters
    ----------
    atom_number: int
        atom number
    fragment: Fragment
        fragment
    brics_fragments: list(Fragment)
        list of BRICS fragments of the ligand as Fragment objects

    Returns
    -------
    tuple(int, str, int/str)
        atom's number, subpocket and BRICS environment
    """

    # get the atom's subpocket (= the fragment's subpocket)
    subpocket = fragment.subpocket.name

    # get environment type of the brics fragment that the current atom belongs to
    environment = next(
        brics_fragment.environment
        for brics_fragment
        in brics_fragments
        if atom_number in brics_fragment.atomNumbers
    )

    return atom_number, subpocket, environment


def _get_atom_properties_dummy(atom_number, atom, atom_tuples, fragments):

    """
    Set atom properties for an atom which is a dummy atom:
    each atom's number, subpocket and BRICS environment.

    Difficulty:
    Dummy atoms have been consecutively numbered as additional ligand atoms.
    However, needed is the atom number of the atom, which was replaced by the dummy atom (dummy-replaced atom).
    The dummy-replaced atom is backtracked via its bond to the atom connected to the dummy atom (dummy-connected atom).
    Note, the dummy-connected atom lives in the same fragment as the dummy atom, whereas the dummy-replaced atom lives
    in another, neighboring fragment.

    Parameters
    ----------
    atom_number: int
        atom number
    atom: Chem.rdchem.Atom
        dummy atom
    atom_tuples: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand
        (final bonds, NOT BRICS bonds!)
    fragments: list(Fragment)
        list of fragments that the ligand consists of

    Returns
    -------
    tuple(int, str, int/str)
        atom's number, subpocket and BRICS environment
    """

    if atom.GetSymbol() != '*':
        raise ValueError(f'Input atom must be dummy atom.')

    # get dummy-connected atom (int)
    dummy_connected_atom_number = _get_dummy_connected_atom_from_dummy_atom(atom)

    # get all bond(s) (atom tuple(s)), which include the dummy-connected atom number (list(tuple(int)))
    bonds_involving_dummy_connected_atom = _get_bonds_from_dummy_connected_atom(
        dummy_connected_atom_number,
        atom_tuples
    )

    # get all dummy-replaced atoms, i.e. atoms that were replaced by the dummy atom(s) during fragmentation (list(int))
    dummy_replaced_atom_numbers = _get_dummy_replaced_atoms_from_bonds(
        dummy_connected_atom_number,
        bonds_involving_dummy_connected_atom
    )

    # if multiple dummy atoms at dummy-connecting atom, select the correct dummy-replaced atom number (int)
    dummy_replaced_atom_number = _select_dummy_replaced_atom(dummy_replaced_atom_numbers, atom_number, fragments)

    # get subpocket of the dummy-replaced atom
    subpocket = _get_dummy_replaced_atom_subpocket(dummy_replaced_atom_number, fragments)

    # dummy atoms do not get an environment type assigned
    environment = 'na'

    return dummy_replaced_atom_number, subpocket, environment


def _get_dummy_connected_atom_from_dummy_atom(atom):
    """
    Get dummy-connected atom number, i.e. the atom that is connected to the dummy atom (all in the same fragment).

    Parameters
    ----------
    atom: Chem.rdchem.Atom
        dummy atom

    Returns
    -------
    int
        atom number of dummy-connected atom
    """

    if atom.GetSymbol() != '*':
        raise ValueError(f'Input atom is not dummy atom.')

    # get neighbor atom connected to dummy atoms (dummy-connected atom)
    dummy_connected_atoms = atom.GetNeighbors()

    # this should yield only one dummy-connected atom since BRICS does not cut in rings
    if len(dummy_connected_atoms) == 1:
        dummy_connected_atom = dummy_connected_atoms[0]
    else:
        raise ValueError(f'Unexpected number of dummy-connected atoms: {len(dummy_connected_atoms)}')

    dummy_connected_atom_number = dummy_connected_atom.GetIntProp('atomNumber')

    return dummy_connected_atom_number


def _get_bonds_from_dummy_connected_atom(dummy_connected_atom_number, atom_tuples):
    """
    Get all ligand cleavage bonds involving the input dummy-connected atom.

    Parameters
    ----------
    dummy_connected_atom_number: int
        atom number of atom connected to dummy atom (dummy-connected atom)
    atom_tuples: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand
        (final bonds, NOT BRICS bonds!)

    Returns
    -------
    list(tuple(int))
        list of atom index tuples which involve the neighbor atom
    """

    # select only ligand cleavage bonds involving the dummy-connected atom
    bonds_involving_dummy_connected_atom = [
        atom_tuple
        for atom_tuple
        in atom_tuples
        if dummy_connected_atom_number in atom_tuple
    ]

    return bonds_involving_dummy_connected_atom


def _get_dummy_replaced_atoms_from_bonds(dummy_connected_atom_number, atom_tuples_involving_neighbor):
    """
    Get the number(s) of dummy-replaced atoms, i.e. of atoms that were replaced by dummy atoms during fragmentation.

    Parameters
    ----------
    dummy_connected_atom_number: int
        atom number of atom connected to dummy atom (dummy-connected atom)
    atom_tuples_involving_neighbor: list(tuple(int))
        list of atom index tuples which involve the neighbor atom

    Returns
    -------
    list(int)
        List of atom numbers that are connected to neighbor atom in full ligand.
    """

    dummy_replaced_atom_numbers = []

    for atom_tuple_involving_neighbor in atom_tuples_involving_neighbor:

        # Select dummy-replaced atom from (dummy-replaced atom, dummy-connected atom) tuple/pair
        dummy_replaced_atom_number = [
            atom_number
            for atom_number
            in atom_tuple_involving_neighbor
            if atom_number != dummy_connected_atom_number
        ]

        if len(dummy_replaced_atom_number) != 1:
            raise ValueError(f'Unexpected number of dummy-replaced atoms: {len(dummy_replaced_atom_number)}')

        dummy_replaced_atom_numbers.extend(
           dummy_replaced_atom_number
        )

    return dummy_replaced_atom_numbers


def _select_dummy_replaced_atom(dummy_replaced_atom_numbers, dummy_atom_number, fragments):
    """
    Select dummy-replaced atom number:
    If only one dummy-replaced atom, select this dummy-replaced atom (this is almost always the case).
    If multiple dummy-replaced atoms are available (this is rarely the case when one atom is connected to multiple
    dummy atoms),select the correct dummy-replaced atom based on matching coordinates of the dummy-replaced atom and
    the dummy atom.

    Parameters
    ----------
    dummy_replaced_atom_numbers: list(int)
        dummy-replaced atom numbers (number retrieved from atom which was replaced by dummy atom during fragmentation)
    dummy_atom_number: int
        dummy atom number (initial number cast during fragmentation)
    fragments: list(Fragment)
        list of fragments that the ligand consists of

    Returns
    -------
    int
        dummy-replaced atom number whose coordinates match the dummy atom's coordinates
    """

    # in almost all cases one dummy atom per dummy-connected atom
    if len(dummy_replaced_atom_numbers) == 1:
        selected_dummy_replaced_atom_number = dummy_replaced_atom_numbers[0]

    # in very rare cases more than one dummy atom per dummy-connected atom
    else:

        selected_dummy_replaced_atom_number = []  # Should have 1 element at the end of iteration

        dummy_atom_position = _get_atom_position_from_atom_number(dummy_atom_number, fragments)

        for dummy_replaced_atom_number in dummy_replaced_atom_numbers:
            dummy_replaced_atom_position = _get_atom_position_from_atom_number(dummy_replaced_atom_number, fragments)

            if np.isclose(dummy_replaced_atom_position.Length(), dummy_atom_position.Length(), rtol=1e-04):
                selected_dummy_replaced_atom_number.append(dummy_replaced_atom_number)

        if len(selected_dummy_replaced_atom_number) != 1:
            raise ValueError(f'Unexpected number of selected dummy-replaced atoms: '
                             f'{len(selected_dummy_replaced_atom_number)}')

        selected_dummy_replaced_atom_number = selected_dummy_replaced_atom_number[0]

    return selected_dummy_replaced_atom_number


def _get_atom_position_from_atom_number(atom_number_of_interest, fragments):

    """
    Get atom position based on the atom number (iterate over fragments to fetch correct atom).

    Parameters
    ----------
    atom_number_of_interest: int
        atom number of interest
    fragments: list(Fragment)
        list of fragments that the ligand consists of

    Returns
    -------
    rdkit.Geometry.rdGeometry.Point3D
        3D coordinates of atom of interest
    """

    atom_position = []  # Should have 1 element at the end of iteration

    for fragment in fragments:

        for i, atom_number in enumerate(fragment.atomNumbers):

            if atom_number == atom_number_of_interest:
                atom_position.append(fragment.mol.GetConformer().GetAtomPosition(i))

    if len(atom_position) != 1:
        raise ValueError(f'Unexpected number of atom positions: {len(atom_position)}')

    return atom_position[0]


def _get_dummy_replaced_atom_subpocket(dummy_replaced_atom_number, fragments):

    """
    Get the dummy atom's subpocket; search in all fragments for the dummy atom number.
    Parameters
    ----------
    dummy_replaced_atom_number: int
        dummy-replaced atom number
    fragments: list(Fragment)
        list of fragments that the ligand consists of

    Returns
    -------
    str
        subpocket name
    """

    dummy_replaced_atom_subpockets = [
        fragment.subpocket
        for fragment
        in fragments
        if dummy_replaced_atom_number in fragment.atomNumbers
    ]

    # this should yield only one subpocket atom since each atom number occurs only once in a ligand
    if len(dummy_replaced_atom_subpockets) == 1:
        dummy_replaced_atom_subpocket = dummy_replaced_atom_subpockets[0]
    else:
        raise ValueError(f'Unexpected number of subpockets: {len(dummy_replaced_atom_subpockets)}')

    return dummy_replaced_atom_subpocket.name
