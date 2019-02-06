from rdkit import Chem
import numpy as np
import pandas as pd
# from rdkit.Chem import Draw
import time

from pocketIdentification import getNearestResidues, getRegion, getSubpocket
from functions import loadAtomInfoFromMol2


# load ligand and binding pocket to rdkit molecules
ligand = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/ligand.mol2', removeHs=False)
pocket = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/pocket.mol2', removeHs=False)
# get molecule conformers
ligandConf = ligand.GetConformer()
pocketConf = pocket.GetConformer()

# read atom information from binding pocket mol2 file (necessary for residue information)
pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/pocket.mol2')

start1 = time.time()

for atom in range(ligand.GetNumAtoms()):
    # get nearest pocket residues
    nearestResidues = getNearestResidues(atom, ligandConf, pocketConf, pocketMol2)
    # get corresponding pocket regions
    regions = [getRegion(res) for res in nearestResidues]


end1 = time.time()
print("loop:", end1 - start1)
# print("Number of distances:", count)

# print(pocketMol2)
