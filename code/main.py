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
for fragIdx, fragment in enumerate(BRICSFragments):

    fragmentAtoms = BRICSFragmentsAtoms[fragIdx]

    # setting atom numbers of fragment as property
    fragment.SetProp('atomNumbers', str(fragmentAtoms))

    subpockets = []

    # get subpocket for each fragment atom
    for i, a in enumerate(fragmentAtoms):

        # the same atom once in the ligand and once in the fragment
        ligandAtom = ligand.GetAtomWithIdx(a)
        fragmentAtom = fragment.GetAtomWithIdx(i)

        # set atom number within the entire molecule as property of the fragment atom
        fragmentAtom.SetProp('atomNumber', str(a))
        # get kinase subpocket that the given atom lies in
        subpocket = getSubpocketFromAtom(a, ligandConf, pocketConf, pocketMol2)
        # ligandAtom.SetProp('subpocket', subpocket)
        subpockets.append(subpocket)

    # get most common subpocket of all atoms in the fragment
    subpocket = mostCommon(subpockets)
    fragment.SetProp('subpocket', subpocket)

    print(fragmentAtoms)
    for bond in fragment.GetBonds():

        endAtom = bond.GetEndAtom()
        beginAtom = bond.GetBeginAtom()

        # if this bond was cleaved
        if endAtom.GetSymbol() == '*':
            print(bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol())

            # find bond in ligand
            beginAtomIdx = beginAtom.GetProp('atomNumber')
            #endAtomIdx = endAtom.GetProp('atomNumber')
            #ligandBond = ligand.GetBondBetweenAtoms(beginAtomIdx, )

       # try:
        #    print(bond.GetBeginAtom().GetProp('atomNumber'), bond.GetEndAtom().GetProp('atomNumber'))
       # except:
       #     print(bond.GetBeginAtom().GetProp('atomNumber'), bond.GetEndAtom().GetSymbol())



    # NEXT:
    # - assign subpockets to fragments
    #       - Do I need fragment + fragmentAtoms?
    #       - fragments as class (extending rdkit molecule class) or just setting properties?
    # - check validity of neighboring fragments

    #       - how to know which fragments are next to each other?
    #       --- store information when finding fragments already?

    # - fragment at certain positions
    #       - which fragmentation function to use?? or recombine fragments?

end1 = time.time()
print("loop:", end1 - start1)
