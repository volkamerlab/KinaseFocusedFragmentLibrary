import numpy as np
import sys

from functions import remove_duplicates, calc_3d_dist, get_ca_atom
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
    String
        Name of the subpocket

    """

    smallest_distance = sys.maxsize  # set smallest distance as max integer value
    nearest_subpocket = Subpocket('noSubpocket')
    for subpocket in subpockets:
        distance = calc_3d_dist(pos, subpocket.center)
        if distance < smallest_distance:
            nearest_subpocket = subpocket
            smallest_distance = distance

    return nearest_subpocket.name


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
            neighboring_fragment = [fragment for fragment in fragments if bond[1] in fragment.atomNumbers][0]
            neighboring_fragments.append(neighboring_fragment)
        elif bond[1] in fragment.atomNumbers:
            neighboring_fragment = [fragment for fragment in fragments if bond[0] in fragment.atomNumbers][0]
            neighboring_fragments.append(neighboring_fragment)
    return neighboring_fragments


def fix_small_fragments(fragments, bonds):

    """
    Fix subpocket assignment of small BRICS fragments such that the final result does not have single small fragments in one subpocket (in place)

    Parameters
    ----------
    fragments: list(Fragment)
        all BRICS fragments of the ligand
    bonds: list(tuple(int))
        list of atom index tuples, where each tuple represents a bond between two atoms in the ligand

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
            if BRICSFragment.mol.GetNumHeavyAtoms() <= 3:
                # find neighboring fragments and fragments in same subpocket
                neighboring_fragments = find_neighboring_fragments(BRICSFragment, fragments, bonds)
                fragments_in_same_subpocket = [f for f in neighboring_fragments if f.subpocket == BRICSFragment.subpocket]
                fragment_sizes = [f.mol.GetNumHeavyAtoms() for f in neighboring_fragments]

                # if small fragment is not yet connected to another fragment
                if not fragments_in_same_subpocket:
                    # connect fragment to largest neighboring fragment (this will also fix single terminal fragments)
                    BRICSFragment.subpocket = neighboring_fragments[int(np.argmax(fragment_sizes))].subpocket
                    num_fixed_fragments += 1

                # if small fragment is already connected to other fragments
                else:
                    subpocket_size = BRICSFragment.mol.GetNumHeavyAtoms() + sum([f.mol.GetNumHeavyAtoms() for f in fragments_in_same_subpocket])
                    # if those fragments build up a large enough fragment, do nothing
                    if subpocket_size > 3:
                        continue
                    # else check further neighboring fragments
                    # (1 round is enough because that would make 3 fragments which should always have a combined size of > 3)
                    else:
                        # find more fragments in the same subpocket
                        for fragment_2 in fragments_in_same_subpocket:
                            fragments_in_same_subpocket_2 = [f for f in find_neighboring_fragments(fragment_2, fragments, bonds)
                                                             if f.subpocket == BRICSFragment.subpocket]
                            # if fragment has neighbors in this subpocket other than BRICSFragment
                            if len(fragments_in_same_subpocket_2) > 1:
                                subpocket_size += (sum([f.mol.GetNumHeavyAtoms() for f in fragments_in_same_subpocket_2])
                                                   - fragment_2.mol.GetNumHeavyAtoms())

                        # if combined fragments in this subpocket are large enough, do nothing
                        if subpocket_size > 3:
                            continue
                        else:
                            # if this is a terminal fragment, do nothing
                            if len(neighboring_fragments) == 1:
                                continue
                            else:
                                # else connect fragment to largest neighboring fragment in other pocket
                                fragments_in_other_subpocket = [f for f in neighboring_fragments if f.subpocket != BRICSFragment.subpocket]
                                fragment_sizes = [f.mol.GetNumHeavyAtoms() for f in fragments_in_other_subpocket]
                                BRICSFragment.subpocket = fragments_in_other_subpocket[int(np.argmax(fragment_sizes))].subpocket
                                num_fixed_fragments += 1
    return None


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


# ================================== OBSOLETE ==========================================

# function that checks validity of neighboring fragments
def checkSubpockets(sp1, sp2):

    subpockets = [sp1, sp2]

    if sp1 == sp2:
        return True
    elif "AP" in subpockets:
        if "FP" in subpockets or "SE" in subpockets or "GA" in subpockets:
            return True
    elif "GA" in subpockets:
        if "FP" in subpockets or "AP" in subpockets or "BP" in subpockets or "B1" in subpockets or "B2" in subpockets:
            return True
    elif "FP" in subpockets:
        if "AP" in subpockets or "GA" in subpockets or "SE" in subpockets:
            return True
    elif "B2" in subpockets:
        if "B1" in subpockets:
            return True
    else:
        return False


# given an atom number of the ligand, get the subpocket that atom lies in
# subpockets: AP, FP, SE, GA, BP
def getSubpocketFromAtom(ligandAtom, ligandConf, subpockets):

    pos = ligandConf.GetAtomPosition(ligandAtom)
    return get_subpocket_from_pos(pos, subpockets)


# given an atom (atom number) of a ligand, find the three nearest protein residues
# pocketMol2: mol2 string of the binding pocket atoms including residue information
#             (all other information in the mol2 string has to be removed)
# ligandAtom: number of the atom of interest in the ligand
def getNearestResidues(ligandAtom, ligandConf, pocketConf, residues):

    lenPocket = pocketConf.GetNumAtoms()
    distances = np.zeros(lenPocket)

    pos1 = ligandConf.GetAtomPosition(ligandAtom)

    # calculate distances from ligand atom to all pocket atoms
    for pocketAtom in range(lenPocket):
        pos2 = pocketConf.GetAtomPosition(pocketAtom)
        distances[pocketAtom] = calc_3d_dist(pos1, pos2)

    # sort pocket atoms by distance
    nearestAtoms = distances.argsort()
    # nearest atoms and corresponding pocket residue numbers
    return remove_duplicates([residues[nearestAtoms[i]] for i in range(lenPocket)])[:3]


# given a list of binding pocket regions, return subpocket
# subpockets: SE, AP, FP, GA, BP,
# SUBPOCKET DEFINITIONS ARE NOT CORRECT!
def getSubpocketFromRegions(regions):

    if 'hinge' in regions[:1]:
        return 'AP'
    elif 'GK' and 'K17' in regions[:1]:
        return 'GA'
    elif 'linker' in regions[:1]:
        return 'SE'
    elif 'linker' and 'DFG' in regions:
        return 'FP'
    elif 'alphaC' in regions:
        return 'BP'
    else:
        return 'other'


# given an atom number of the ligand, get the subpocket that atom lies in
# subpockets: AP, FP, SE, GA, BP
def getSubpocketFromAtomDistances(ligandAtom, ligandConf, pocketConf, residues):

    # get nearest pocket residues
    nearestResidues = getNearestResidues(ligandAtom, ligandConf, pocketConf, residues)

    # get corresponding pocket regions
    regions = [get_region(res) for res in nearestResidues]

    return getSubpocketFromRegions(regions)