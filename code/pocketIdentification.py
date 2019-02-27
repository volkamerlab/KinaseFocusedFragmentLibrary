import numpy as np
import sys

from functions import removeDuplicates, calculate3DDistance


# given an atom (atom number) of a ligand, find the three nearest protein residues
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
    elif res == 17:
        return 'K17'
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


# given a list of binding pocket regions, return subpocket
# subpockets: SE, AP, FP, GA, BP,
def getSubpocketFromRegions(regions):

    if 'hinge' in regions[:1]:
        return 'AP'
    elif 'GK' and 'K17' in regions:
        return 'GA'
    elif 'linker' and 'DFG' in regions:
        return 'FP'
    elif 'alphaC' in regions:
        return 'BP'
    else:
        return 'other'


# given an atom number of the ligand, get the subpocket that atom lies in
# subpockets: AP, FP, SE, GA, BP
def getSubpocketFromAtom(ligandAtom, ligandConf, pocketConf, pocketMol2):

    # get nearest pocket residues
    nearestResidues = getNearestResidues(ligandAtom, ligandConf, pocketConf, pocketMol2)

    # get corresponding pocket regions
    regions = [getRegion(res) for res in nearestResidues]

    return getSubpocketFromRegions(regions)


# function that checks validity of neighboring fragments
def checkSubpockets(sb1, sb2):

    subpockets = [sb1, sb2]

    if sb1 == sb2:
        return True
    elif "AP" in subpockets:
        if "FP" or "SE" or "GA" in subpockets:
            return True
        else:
            return False
    elif "GA" in subpockets:
        if "FP" or "AP" or "BP" in subpockets:
            return True
        else:
            return False
    elif "other" in subpockets:
        return True
    else:
        return False
