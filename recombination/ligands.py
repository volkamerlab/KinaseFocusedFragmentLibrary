from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.PropertyMol import PropertyMol
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)

from pathlib import Path
import pickle
import sys
import time
sys.path.append("../fragmentation/")

from pickle_loader import pickle_loader
from construct_ligand import construct_ligand

start = time.time()

# ============================= LIGAND CONSTRUCTION ============================================

in_path = Path('meta_library.pickle')
pickle_in = in_path.open('rb')
out_file = Path('../CombinatorialLibrary/combinatorial_library.pickle').open('wb')

ligand_smiles = set()
# iterate over ligands
for meta in pickle_loader(pickle_in):

    ligand = construct_ligand(meta, ligand_smiles=ligand_smiles)

    pickle.dump(PropertyMol(ligand), out_file)
    ligand_smiles.add(Chem.MolToSmiles(ligand))

out_file.close()

runtime = time.time() - start
print('Number of resulting ligands: ', len(ligand_smiles))
print('Time: ', runtime)

with open('ligand_smiles.pickle', 'wb') as pickle_out:
    for smiles in ligand_smiles:
        pickle.dump(smiles, pickle_out)



