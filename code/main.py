from rdkit import Chem
from rdkit.Chem import Draw
import time
import sys

from pocketIdentification import getSubpocketFromAtom, checkSubpockets
from functions import loadAtomInfoFromMol2, mostCommon
from fragmentation import FindBRICSFragments, GetFragmentsFromAtomTuples
from classes import Fragment

# load ligand and binding pocket to rdkit molecules
ligand = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/ligand.mol2', removeHs=False)
pocket = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2', removeHs=False)

# get molecule conformers
ligandConf = ligand.GetConformer()
pocketConf = pocket.GetConformer()

lenLigand = ligand.GetNumAtoms()

# read atom information from binding pocket mol2 file (necessary for residue information)
pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')

start = time.time()

# get subpocket for each ligand atom
for a, atom in enumerate(ligand.GetAtoms()):

    subpocket = getSubpocketFromAtom(a, ligandConf, pocketConf, pocketMol2)
    atom.SetProp('subpocket', subpocket)

end = time.time()
print("Subpocket identification:", end - start)
start = time.time()


# find BRICS fragments and bonds
BRICSFragmentsAtoms, BRICSBonds = FindBRICSFragments(ligand)

# list to store the bonds where we will cleave
bonds = []
# BRICS fragments as Fragment objects
BRICSFragments = [Fragment(atomNumbers=BRICSFragmentsAtoms[f]) for f in range(len(BRICSFragmentsAtoms))]

# iterate over BRICS bonds
for beginAtom, endAtom in BRICSBonds:

    # find corresponding fragments
    firstFragment = [fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers][0]
    secondFragment = [fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers][0]

    # add subpocket to fragment objects (if not yet defined for this fragment)
    if firstFragment.subpocket is None:
        firstSubpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in firstFragment.atomNumbers])
        firstFragment.subpocket = firstSubpocket
    else:
        firstSubpocket = firstFragment.subpocket
    if secondFragment.subpocket is None:
        secondSubpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in secondFragment.atomNumbers])
        secondFragment.subpocket = secondSubpocket
    else:
        secondSubpocket = secondFragment.subpocket

    # check validity of subpockets
    if not checkSubpockets(firstSubpocket, secondSubpocket):

        print("ERROR: Subpockets "+firstSubpocket+" and "+secondSubpocket+" can not be connected. "
                                                                          "Molecule is skipped.")
        # skip this molecule
        sys.exit()  # change this line when we work with more molecules -> continue with for loop

    # if subpockets of the 2 fragments differ (and they are valid)
    if firstSubpocket != secondSubpocket:

        # store this bond as a bond where we will cleave
        bonds.append((beginAtom, endAtom))


# actual fragmentation
fragments = GetFragmentsFromAtomTuples(bonds, BRICSFragments, ligand)

img = Draw.MolsToGridImage([fragment.mol for fragment in fragments],
                           legends=[fragment.subpocket for fragment in fragments])
img.save('test/3w2s_subpocket-fragmentation.png')


end = time.time()
print("Fragmentation:", end - start)


# TO DO:
# - store bond information (neighboring subpocket, rule?)
#
# - implement correct subpocket definition (e.g. subpocket centers)
#
# - iteration over multiple molecules
