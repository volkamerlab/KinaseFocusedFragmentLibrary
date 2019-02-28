import numpy as np


# calculate 3D distance between two atoms
def calculate3DDistance(conformer1, conformer2, atom1, atom2):

    pos1 = conformer1.GetAtomPosition(atom1)
    pos2 = conformer2.GetAtomPosition(atom2)
    return np.linalg.norm(pos1 - pos2)


# load atom information block from a mol2 file into a list of lists of strings
def loadAtomInfoFromMol2(file):

    atomInfo = []
    with open(file) as f:
        line = f.readline()
        while not line.startswith('@<TRIPOS>ATOM'):
            line = f.readline()
        line = f.readline()
        while not line.startswith('@<TRIPOS>'):
            atomInfo.append(line.split())
            line = f.readline()

    return atomInfo


# Eva's function
def get_cavity_center_of_mass_ca(res_list):
    """
    For the given residues compute the coordinates of their
    center of mass with respect to their C alpha atoms.

    Input: list of pdb residue entities
    Output: coordinates of CoM (array 1x3)

    """
    # print 'calculating coordinates of current reference location'
    center = np.zeros(3, float)
    tmp = []
    for res in res_list:
        for atom in res:
            # get ca atom
            if atom.get_name() == 'CA':
                tmp.append(atom.get_coord())
                center += atom.get_coord()
    return sum(tmp) / len(res_list)


# remove duplicates from array/list lst
def removeDuplicates(lst):
    return list(dict.fromkeys(lst))


# find most common element in a list
def mostCommon(lst):
    return max(set(lst), key=lst.count)
