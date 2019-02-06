import numpy as np
import time

from functions import removeDuplicates


# calculate 3D distance between two atoms
def calculate3DDistance(conformer1, conformer2, atom1, atom2):

    start = time.time()
    pos1 = conformer1.GetAtomPosition(atom1)
    pos2 = conformer2.GetAtomPosition(atom2)
    end = time.time()
    # print("conformer:", end - start)
    start = time.time()
    result = np.linalg.norm(pos1 - pos2)
    end = time.time()
    # print("distance:", end - start)

    return result


# given an atom of the ligand, find the three nearest protein residues
# pocketMol2: mol2 string of the binding pocket atoms including residue information
#             (all other information in the mol2 string has to be removed)
# ligandAtom: number of the atom of interest in the ligand
def getNearestResidues(ligandAtom, ligandConf, pocketConf, pocketMol2):

    lenPocket = pocketConf.GetNumAtoms()
    distances = np.zeros(lenPocket)
    # calculate distances from ligand atom to all pocket atoms
    for pocketAtom in range(lenPocket):
        distances[pocketAtom] = calculate3DDistance(ligandConf, pocketConf, ligandAtom, pocketAtom)

    # sort pocket atoms by distance
    nearestAtoms = distances.argsort()
    # nearest atoms and corresponding pocket residue numbers
    return removeDuplicates([int(pocketMol2[nearestAtoms[i]-1][6]) for i in range(lenPocket)])[:3]


# given a residue number within the binding pocket (KLIFS numbering):
# returns the corresponding region of the binding pocket (KLIFS definition) as a string
def getRegion(res):

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


# given a list of kinase binding pocket regions, get the corresponding subpocket
# subpockets: AP, FP, SE, GA, BP
# ----------------- TO DO ? ---------------------
def getSubpocket(regions):

    if 'hinge' and 'linker' in regions:
        return 'AP'
    elif 'GK' in regions:
        return 'GA'
    else:
        return 'other'
# -----------------------------------------------
