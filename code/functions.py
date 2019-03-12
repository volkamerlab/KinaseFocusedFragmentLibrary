import numpy as np


def calculate3DDistance(pos1, pos2):

    return np.linalg.norm(pos1 - pos2)


# get CA atom object of a residue number in a pocket
def getCaAtom(res, pocketMol2, pocket):

    # if res not in pocketMol2.res_id.values:
    #    print('ERROR: Important residue is missing in structure. Molecule is skipped.')
    #    return None
    pocketMol2Res = pocketMol2[pocketMol2.res_id == res]
    CaAtom = pocketMol2Res[pocketMol2Res.atom_name == 'CA'].index.values
    if len(CaAtom) == 0:
        print('ERROR: Important residue is missing in structure. Molecule is skipped.')
        return None
    CaAtomId = int(CaAtom[0])
    CaAtom = pocket.GetAtomWithIdx(CaAtomId)
    return CaAtom


# remove duplicates from array/list lst
def removeDuplicates(lst):
    return list(dict.fromkeys(lst))


# find most common element in a list
def mostCommon(lst):
    return max(set(lst), key=lst.count)


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
