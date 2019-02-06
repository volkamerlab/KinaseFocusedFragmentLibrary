from rdkit import Chem
# from rdkit.Chem import BRICS
from rdkit.Chem import Draw

from BRICS import BRICSDecompose

ligand = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/ligand.mol2', removeHs=False)

# res = [Chem.MolFromSmiles(x) for x in list(BRICS.BRICSDecompose(ligand, minFragmentSize=5))]
res = [Chem.MolFromSmiles(x) for x in list(BRICSDecompose(ligand, minFragmentSize=1, silent=False))]

img = Draw.MolsToGridImage(res)
img.save('test/3w2s_ligand_brics.png')
