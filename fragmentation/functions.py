import numpy as np


# get CA atom object of a residue number in a pocket
def get_ca_atom(res, pocket_mol2, pocket):

    """
    Get the C alpha atom of a given residue in the protein/binding pocket

    Parameters
    ----------
    res: int
        residue ID
    pocket_mol2: Pandas DataFrame
        atom block from mol2 file representing the pocket (read using PandasMol2().read_mol2())
    pocket: RDKit Mol object
        protein/binding pocket molecule

    Returns
    -------
    ca_atom: RDKit Atom object or None
        C alpha atom, None if atom or residue was not found

    """

    pocket_mol2_res = pocket_mol2[pocket_mol2.res_id == res]
    ca_atom = pocket_mol2_res[pocket_mol2_res.atom_name == 'CA'].index.values
    if len(ca_atom) == 0:
        # if this residue/atom is missing in the structure, return None
        return None
    ca_atom_id = int(ca_atom[0])
    ca_atom = pocket.GetAtomWithIdx(ca_atom_id)
    return ca_atom


# calculates the distance between two points in 3D
def calc_3d_dist(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


# remove duplicates from array/list lst
def remove_duplicates(lst):
    return list(dict.fromkeys(lst))


# find most common element in a list
def most_common(lst):
    return max(set(lst), key=lst.count)
