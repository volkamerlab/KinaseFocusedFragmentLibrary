# Construct all ligands from Combination objects in file 'meta_library.pickle'

from rdkit import Chem
from rdkit.Chem.PropertyMol import PropertyMol
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)

from pathlib import Path
import pickle
import time
import sys
sys.path.append('../')

from pickle_loader import pickle_loader
from construct_ligand import construct_ligand

start = time.time()

# ============================= LIGAND CONSTRUCTION ============================================

path_to_library = Path('../results')
in_paths = list(path_to_library.glob('*.pickle'))

out_path = '../../CombinatorialLibrary/'
out_file = Path(out_path+'combinatorial_library.pickle').open('wb')

ligand_smiles = set()

# iterate over ligands
for in_path in in_paths:

    with open(in_path, 'rb') as pickle_in:
        for meta in pickle_loader(pickle_in):

            ligand = construct_ligand(meta, ligand_smiles=ligand_smiles)
            if not ligand:
                continue
            property_ligand = PropertyMol(ligand)
            pickle.dump(property_ligand, out_file)
            ligand_smiles.add(Chem.MolToSmiles(ligand))

out_file.close()

runtime = time.time() - start
print('Number of resulting ligands: ', len(ligand_smiles))
print('Time: ', runtime)

with open(out_path+'ligand_smiles.pickle', 'wb') as pickle_out:
    for smiles in ligand_smiles:
        pickle.dump(smiles, pickle_out)
