from rdkit import Chem
# from rdkit.Chem import Draw
import time

from pocketIdentification import getSubpocketFromAtom
from functions import loadAtomInfoFromMol2, mostCommon
from fragmentation import FindBRICSFragments


# load ligand and binding pocket to rdkit molecules
ligand = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/ligand.mol2', removeHs=False)
pocket = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2', removeHs=False)

# get molecule conformers
ligandConf = ligand.GetConformer()
pocketConf = pocket.GetConformer()

lenLigand = ligand.GetNumAtoms()

# read atom information from binding pocket mol2 file (necessary for residue information)
pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')

start1 = time.time()

# find BRICS fragments
BRICSFragmentsAtoms, BRICSFragments = FindBRICSFragments(ligand, asMols=True)

# get subpocket for each BRICS fragment
for i, fragment in enumerate(BRICSFragments):

    fragmentAtoms = BRICSFragmentsAtoms[i]

    subpockets = []

    # get subpocket for each fragment atom
    for a in fragmentAtoms:
        atom = ligand.GetAtomWithIdx(a)
        # get kinase subpocket that the given atom lies in
        subpocket = getSubpocketFromAtom(a, ligandConf, pocketConf, pocketMol2)
        atom.SetProp('subpocket', subpocket)
        subpockets.append(subpocket)

    # get most common subpocket of all atoms in the fragment
    subpocket = mostCommon(subpockets)
    print(subpocket)

    # NEXT:
    # - assign subpockets to fragments
    #       - Do I need fragment + fragmentAtoms?
    #       - fragments as class (extending rdkit molecule class) or just setting properties
    # - check validity of neighboring fragments
    #       - how to know which fragments are next to each other?
    # - fragment at certain positions
    #       - which fragmentation function to use?? or recombine fragments?

end1 = time.time()
print("loop:", end1 - start1)
