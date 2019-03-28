import numpy as np
import sys

from functions import removeDuplicates, calculate3DDistance
from classes import Subpocket


# given an atom number of the ligand, get the subpocket that atom lies in
# subpockets: AP, FP, SE, GA, BP
def getSubpocketFromAtom(ligandAtom, ligandConf, subpockets):

    pos = ligandConf.GetAtomPosition(ligandAtom)
    return getSubpocketFromPos(pos, subpockets)


def getSubpocketFromPos(pos, subpockets):

    smallestDistance = sys.maxsize  # set smallest distance as max integer value
    nearestSubpocket = Subpocket('noSubpocket')
    for subpocket in subpockets:
        distance = calculate3DDistance(pos, subpocket.center)
        if distance < smallestDistance:
            nearestSubpocket = subpocket
            smallestDistance = distance

    return nearestSubpocket.name


def findNeighboringFragments(fragment, fragments, bonds):
    neighboringFragments = []
    for bond in bonds:
        if bond[0] in fragment.atomNumbers:
            neighboringFragment = [fragment for fragment in fragments if bond[1] in fragment.atomNumbers][0]
            neighboringFragments.append(neighboringFragment)
        elif bond[1] in fragment.atomNumbers:
            neighboringFragment = [fragment for fragment in fragments if bond[0] in fragment.atomNumbers][0]
            neighboringFragments.append(neighboringFragment)
    return neighboringFragments


# fix subpockets of small BRICS fragments such that the final result does not have single small fragments
def fixSmallFragments(BRICSFragments, BRICSBonds):

    # repeat until all small fragments are fixed
    numFixedFragments = 1
    while numFixedFragments > 0:
        numFixedFragments = 0

        # iterate over BRICS fragments, which now all have a subpocket assigned
        for BRICSFragment in BRICSFragments:

            # small fragments
            if BRICSFragment.mol.GetNumHeavyAtoms() <= 3:
                # find neighboring fragments and fragments in same subpocket
                neighboringFragments = findNeighboringFragments(BRICSFragment, BRICSFragments, BRICSBonds)
                fragmentsInSameSubpocket = [f for f in neighboringFragments if f.subpocket == BRICSFragment.subpocket]
                fragmentSizes = [f.mol.GetNumHeavyAtoms() for f in neighboringFragments]

                # if small fragment is not yet connected to another fragment
                if not fragmentsInSameSubpocket:
                    # connect fragment to largest neighboring fragment (this will also fix single terminal fragments)
                    BRICSFragment.subpocket = neighboringFragments[int(np.argmax(fragmentSizes))].subpocket
                    numFixedFragments += 1

                # if small fragment is already connected to other fragments
                else:
                    subpocketSize = BRICSFragment.mol.GetNumHeavyAtoms() + sum([f.mol.GetNumHeavyAtoms() for f in fragmentsInSameSubpocket])
                    # if those fragments build up a large enough fragment, do nothing
                    if subpocketSize > 3:
                        continue
                    # else check further neighboring fragments
                    # (1 round is enough because that would make 3 fragments which should always have a combined size of > 3)
                    else:
                        # find more fragments in the same subpocket
                        for fragment in fragmentsInSameSubpocket:
                            fragmentsInSameSubpocket2 = [f for f in findNeighboringFragments(fragment, BRICSFragments, BRICSBonds)
                                                         if f.subpocket == BRICSFragment.subpocket]
                            # if fragment has neighbors in this subpocket other than BRICSFragment
                            if len(fragmentsInSameSubpocket2) > 1:
                                subpocketSize += (sum([f.mol.GetNumHeavyAtoms() for f in fragmentsInSameSubpocket2])-fragment.mol.GetNumHeavyAtoms())

                        # if combined fragments in this subpocket are large enough, do nothing
                        if subpocketSize > 3:
                            continue
                        else:
                            # if this is a terminal fragment, do nothing
                            if len(neighboringFragments) == 1:
                                continue
                            else:
                                # else connect fragment to largest neighboring fragment in other pocket
                                fragmentsInOtherSubpocket = [f for f in neighboringFragments if f.subpocket != BRICSFragment.subpocket]
                                fragmentSizes = [f.mol.GetNumHeavyAtoms() for f in fragmentsInOtherSubpocket]
                                BRICSFragment.subpocket = fragmentsInOtherSubpocket[int(np.argmax(fragmentSizes))].subpocket
                                numFixedFragments += 1
    return None


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


# get geometric center of atoms (list of atom objects) in mol
def getGeometricCenter(atoms, molConf):

    center = np.zeros(3, float)
    for atom in atoms:
        pos = molConf.GetAtomPosition(atom.GetIdx())
        center += pos
    return center / len(atoms)


# ================================== OLD APPROACH ==========================================

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
        distances[pocketAtom] = calculate3DDistance(pos1, pos2)

    # sort pocket atoms by distance
    nearestAtoms = distances.argsort()
    # nearest atoms and corresponding pocket residue numbers
    return removeDuplicates([residues[nearestAtoms[i]] for i in range(lenPocket)])[:3]


# given a residue number within the binding pocket (KLIFS numbering):
# returns the corresponding region of the binding pocket (KLIFS definition) as a string
# USED IN PYMOL SCRIPT!
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
    regions = [getRegion(res) for res in nearestResidues]

    return getSubpocketFromRegions(regions)
