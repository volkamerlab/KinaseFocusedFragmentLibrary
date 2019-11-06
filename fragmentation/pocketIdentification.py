import numpy as np
import sys

from functions import calc_3d_dist, get_ca_atom
from classes import Subpocket


# given a 3D position, get the subpocket of the nearest subpocket center
# subpockets: AP, FP, SE, GA, B1, B2
def get_subpocket_from_pos(pos, subpockets):

    """
    Get the subpocket of the nearest subpocket center to a given 3D position

    Parameters
    ----------
    pos: list(float)
        3D position in the kinase binding pocket
    subpockets: list(Subpocket)
        list of subpockets as Subpocket objects

    Returns
    -------
    nearest_subpocket: Subpocket
    smallest_distance: float
        distance from pos to nearest_subpocket center

    """

    smallest_distance = sys.maxsize  # set smallest distance as max integer value
    nearest_subpocket = Subpocket('noSubpocket')
    for subpocket in subpockets:
        distance = calc_3d_dist(pos, subpocket.center)
        if distance < smallest_distance:
            nearest_subpocket = subpocket
            smallest_distance = distance

    return nearest_subpocket, smallest_distance


def find_neighboring_fragments(fragment, fragments, bonds):

    """
    Get all fragments connected of a given fragment

    Parameters
    ----------
    fragment: Fragment
        fragment of interest
    fragments: list(Fragment)
        all fragments of the ligand
    bonds: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand

    Returns
    -------
    neighboring_fragments: list(Fragment)
        list of fragments connected to the given fragment

    """

    neighboring_fragments = []
    for bond in bonds:
        if bond[0] in fragment.atomNumbers:
            neighboring_fragment = next(fragment for fragment in fragments if bond[1] in fragment.atomNumbers)
            neighboring_fragments.append(neighboring_fragment)
        elif bond[1] in fragment.atomNumbers:
            neighboring_fragment = next(fragment for fragment in fragments if bond[0] in fragment.atomNumbers)
            neighboring_fragments.append(neighboring_fragment)
    return neighboring_fragments


def fix_small_fragments(fragments, bonds, min_size):

    """
    Fix subpocket assignment of small BRICS fragments such that the final result does not have single small fragments in one subpocket (in place)

    Parameters
    ----------
    fragments: list(Fragment)
        all BRICS fragments of the ligand
    bonds: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand
    min_size: int
        Minimum allowed fragment size

    Returns
    -------
    None

    """

    # repeat until all small fragments are fixed
    num_fixed_fragments = 1
    while num_fixed_fragments > 0:

        num_fixed_fragments = 0

        # iterate over BRICS fragments, which now all have a subpocket assigned
        for BRICSFragment in fragments:

            # small fragments
            if BRICSFragment.mol.GetNumHeavyAtoms() < min_size:

                # find neighboring fragments
                neighboring_fragments = find_neighboring_fragments(BRICSFragment, fragments, bonds)
                fragments_in_same_subpocket = [f for f in neighboring_fragments if f.subpocket == BRICSFragment.subpocket]
                fragments_in_other_subpocket = [f for f in neighboring_fragments if f.subpocket != BRICSFragment.subpocket]

                # if this fragment has only neighbors within the same subpocket, it can not be fixed in this iteration
                if not fragments_in_other_subpocket:
                    continue

                # if small fragment is not yet connected to another fragment
                elif not fragments_in_same_subpocket:
                    fragment_sizes = [f.mol.GetNumHeavyAtoms() for f in neighboring_fragments]
                    # connect fragment to largest neighboring fragment (this will also fix single terminal fragments)
                    BRICSFragment.subpocket = neighboring_fragments[int(np.argmax(fragment_sizes))].subpocket
                    num_fixed_fragments += 1

                # if small fragment is already connected to other fragments
                else:

                    # if this is a terminal fragment, do nothing
                    if len(neighboring_fragments) == 1:
                        continue

                    # all fragments in this subpocket
                    subpocket_fragments = fragments_in_same_subpocket+[BRICSFragment]
                    n_subpocket_fragments = len(subpocket_fragments)

                    subpocket_size = sum([f.mol.GetNumHeavyAtoms() for f in subpocket_fragments])

                    # look for further neighboring fragments in the same subpocket until min size is reached or all fragments have been found
                    n_start = 1
                    while subpocket_size < min_size and n_subpocket_fragments > n_start:

                        n_start = n_subpocket_fragments

                        # find more fragments in the same subpocket
                        for fragment in fragments_in_same_subpocket:
                            fragments_in_same_subpocket_2 = [f for f in find_neighboring_fragments(fragment, fragments, bonds)
                                                             if f.subpocket == BRICSFragment.subpocket]
                            subpocket_fragments = subpocket_fragments + [f for f in fragments_in_same_subpocket_2
                                                                         if f.atomNumbers not in [g.atomNumbers for g in subpocket_fragments]]
                            n_subpocket_fragments = len(subpocket_fragments)
                            # if fragment has neighbors in this subpocket other than BRICSFragment
                            subpocket_size = sum([f.mol.GetNumHeavyAtoms() for f in subpocket_fragments])
                            # go to next fragment if min size is fulfilled
                            if subpocket_size >= min_size:
                                break

                        fragments_in_same_subpocket = fragments_in_same_subpocket_2

                    # if combined fragments in this subpocket are not large enough
                    if subpocket_size < min_size:
                        # else connect fragment to largest neighboring fragment in other pocket
                        fragment_sizes = [f.mol.GetNumHeavyAtoms() for f in fragments_in_other_subpocket]
                        BRICSFragment.subpocket = fragments_in_other_subpocket[int(np.argmax(fragment_sizes))].subpocket
                        num_fixed_fragments += 1


def calc_geo_center(atoms, mol_conf):

    """
    Calculate geometric center of given atoms in a molecule

    Parameters
    ----------
    atoms: list(Atom)
        list of RDKit Atom objects
    mol_conf: RDKit molecule conformer
        conformer of the molecule containing the atoms

    Returns
    -------
    list(float)
        3D position of the geometric center

    """

    center = np.zeros(3, float)
    for atom in atoms:
        pos = mol_conf.GetAtomPosition(atom.GetIdx())
        center += pos
    return center / len(atoms)


def calc_subpocket_center(subpocket, pocket, pocket_mol2, folder):

    """
    Calculate geometric center of a given subpocket (NOT in place, Subpocket object is not changed)

    Parameters
    ----------
    subpocket: Subpocket
        Subpocket object to calculate the center of
    pocket: RDKit Mol object
        binding pocket molecule
    pocket_mol2: Pandas DataFrame
        atom block from mol2 file representing the pocket (read using PandasMol2().read_mol2())
    folder: String
        Structure ID, e.g. '3w2s_altA_chainA'

    Returns
    -------
    list(float)
        3D position of the geometric center

    """

    pocket_conf = pocket.GetConformer()
    ca_atoms = []
    for res in subpocket.residues:
        ca_atom = get_ca_atom(res, pocket_mol2, pocket)
        # if this residue or its C alpha atom is missing
        if not ca_atom:
            # print(res)
            # try neighboring residues
            ca_atoms_nei = []
            for resNei in [res - 1, res + 1]:
                ca_atoms_nei.append(get_ca_atom(resNei, pocket_mol2, pocket))
            if None in ca_atoms_nei:
                # if only one neighboring residue is missing, take the other one
                ca_atom = [atom for atom in ca_atoms_nei if atom]
                if ca_atom:
                    ca_atom = ca_atom[0]
                    ca_atom_pos = pocket_conf.GetAtomPosition(ca_atom.GetIdx())
                # if both neighboring residues are missing
                else:
                    print('ERROR in ' + folder + ':')
                    print('Important residue ' + str(res) + ' is missing in structure. Structure is skipped. \n')
                    return None
            else:
                # if both neighboring residues are present, take the center
                ca_atom_pos = calc_geo_center(ca_atoms_nei, pocket_conf)
        else:
            ca_atom_pos = pocket_conf.GetAtomPosition(ca_atom.GetIdx())

        ca_atoms.append(ca_atom_pos)

    # overwrite subpocket center for current structure
    center = np.zeros(3, float)
    for pos in ca_atoms:
        center += pos

    return center / len(ca_atoms)


# check validity of neighboring fragments
valid_subpocket_connections = [{'SE', 'AP'},
                               {'SE', 'FP'},
                               {'AP', 'FP'},
                               {'AP', 'GA'},
                               {'FP', 'GA'},
                               {'GA', 'B1'},
                               {'GA', 'B2'},
                               {'B1', 'B2'}]


def is_valid_subpocket_connection(sp_1, sp_2):

    """
    Given two subpocket, checks whether fragments assigned to these subpockets are allowed to be connected

    Parameters
    ----------
    sp_1, sp_2: Subpocket
        two subpocket objects

    Returns
    -------
    True if connecting two fragments assigned to these subpockets is allowed,
    False otherwise.

    """

    if {sp_1.name, sp_2.name} in valid_subpocket_connections or 'X' in sp_1.name or 'X' in sp_2.name:
        return True
    else:
        return False


# given a residue number within the binding pocket (KLIFS numbering):
# returns the corresponding region of the binding pocket (KLIFS definition) as a string
# used only in pymol script
def get_region(res):

    """
    Given a residue number within the binding pocket (KLIFS numbering):
    Returns the corresponding region of the binding pocket (KLIFS definition) as a string

    Parameters
    ----------
    res: int
        Residue ID (1-85)

    Returns
    -------
    String
        Kinase binding pocket region

    """

    if 1 <= res <= 3:
        return 'beta1'
    elif 4 <= res <= 9:
        return 'g.l'  # glycine rich loop
    elif 10 <= res <= 13:
        return 'beta2'
    elif 14 <= res <= 19:
        return 'beta3'
    elif 20 <= res <= 30:
        return 'alphaC'
    elif 31 <= res <= 37:
        return 'b.l'
    elif 38 <= res <= 41:
        return 'beta4'
    elif 42 <= res <= 44:
        return 'beta5'
    elif res == 45:
        return 'GK'
    elif 46 <= res <= 48:
        return 'hinge'
    elif 49 <= res <= 52:
        return 'linker'
    elif 53 <= res <= 59:
        return 'alphaD'
    elif 60 <= res <= 64:
        return 'alphaE'
    elif 65 <= res <= 67:
        return 'beta6'
    elif 68 <= res <= 75:
        return 'c.l'
    elif 76 <= res <= 78:
        return 'beta7'
    elif res == 79:
        return 'beta8'
    elif 80 <= res <= 83:
        return 'DFG'
    elif 84 <= res <= 85:
        return 'a.I'
    else:
        sys.exit('ERROR: Given residue number not between 1 and 85!')
