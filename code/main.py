from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import time
import sys

from pocketIdentification import getSubpocketFromAtom, checkSubpockets
from functions import loadAtomInfoFromMol2, mostCommon
from fragmentation import FindBRICSFragments, getFragmentsFromAtomTuples
from classes import Fragment

from biopandas.mol2 import PandasMol2


# load ligand and binding pocket to rdkit molecules
ligand = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/ligand.mol2', removeHs=False)
pocket = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2', removeHs=False)

# get molecule conformers
ligandConf = ligand.GetConformer()
pocketConf = pocket.GetConformer()

lenLigand = ligand.GetNumAtoms()

# read atom information from binding pocket mol2 file (necessary for residue information)
# pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')
pocketMol2 = PandasMol2().read_mol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2',
                                    columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                             5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                             9: ('secondary structure', str)}).df
residues = pocketMol2.res_id.apply(int)

start = time.time()

# calculate subpocket centers
# - for each subpocket:
#       - get residues defining the pocket
#       - get C alpha atoms of these residues
#       - calculate center of mass of these atoms

# get subpocket for each ligand atom
for a, atom in enumerate(ligand.GetAtoms()):

    subpocket = getSubpocketFromAtom(a, ligandConf, pocketConf, residues)
    atom.SetProp('subpocket', subpocket)

end = time.time()
print("Subpocket identification:", end - start)
start = time.time()


# find BRICS fragments and bonds (as atom numbers)
BRICSFragmentsAtoms, BRICSBonds = FindBRICSFragments(ligand)

# list to store the bonds where we will cleave (as atom tuples)
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
fragments = getFragmentsFromAtomTuples(bonds, BRICSFragments, ligand)


for fragment in fragments:
    tmp = AllChem.Compute2DCoords(fragment.mol)

# atom = fragment.mol.GetAtoms()[0]
# print(atom.GetProp('atomNumber'), atom.GetProp('neighboringSubpocket'), atom.GetProp('subpocket'), atom.GetProp('priority'),
#     [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()])

img = Draw.MolsToGridImage([fragment.mol for fragment in fragments],
                           legends=[fragment.subpocket for fragment in fragments],
                           subImgSize=(400, 400))
img.save('test/3w2s_subpocket-fragmentation.png')


end = time.time()
print("Fragmentation:", end - start)


# TO DO:
# - store bond information (BRICS rule)
#
# - implement correct subpocket definition (subpocket centers)
#
# - iteration over multiple molecules
