from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import Draw

from myBRICS import FindBRICSBonds, BreakBRICSBonds, BRICSDecompose

ligand = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/ligand.mol2', removeHs=False)

# -------- BreakBRICSBonds ---------

res = BreakBRICSBonds(ligand)
res = Chem.MolToSmiles(res).split('.')
res = [Chem.MolFromSmiles(x) for x in res]

img = Draw.MolsToGridImage(res)

img.save('test/3w2s_ligand_BreakBRICSBonds.png')


# -------- FindBRICSBonds ---------

# get bonds that BRICS would cleave -- go from here and discard some bonds? -- probably inefficient...
atomTuples = [bond[0] for bond in list(FindBRICSBonds(ligand))]
bonds = [ligand.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atomTuples]
res = Chem.FragmentOnBonds(ligand, bonds)
res = Chem.MolToSmiles(res).split('.')
res = [Chem.MolFromSmiles(x) for x in res]

img = Draw.MolsToGridImage(res)

img.save('test/3w2s_ligand_FindBRICSBonds.png')


# -------- BRICSDecompose ---------

res = [Chem.MolFromSmiles(x) for x in list(BRICSDecompose(ligand, minFragmentSize=5, silent=True))]
print([Chem.MolToSmiles(x) for x in res])

img = Draw.MolsToGridImage(res)

img.save('test/3w2s_ligand_BRICSDecompose.png')
