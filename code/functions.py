import numpy as np
from collections import Counter


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


# remove duplicates from array/list lst
def removeDuplicates(lst):
    return list(dict.fromkeys(lst))


# find most common element in a list
def mostCommon(lst):
    return max(set(lst), key=lst.count)
