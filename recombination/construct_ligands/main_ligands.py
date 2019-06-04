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
from construct_ligand import construct_ligand, read_fragment_library

start = time.time()

data = read_fragment_library(Path('../../FragmentLibrary'))

# ============================= LIGAND CONSTRUCTION ============================================

path_to_library = Path('../results')
in_paths = path_to_library.glob('*.pickle')

out_path = '../../CombinatorialLibrary/combinatorial_library.pickle'
with open(out_path, 'wb') as out_file:

    ligand_smiles = set()

    # iterate over ligands
    for in_path in in_paths:

        with open(in_path, 'rb') as pickle_in:
            for meta in pickle_loader(pickle_in):

                ligand = construct_ligand(meta, data, ligand_smiles=ligand_smiles)
                if not ligand:
                    continue
                property_ligand = PropertyMol(ligand)
                pickle.dump(property_ligand, out_file)
                ligand_smiles.add(Chem.MolToSmiles(ligand))

runtime = time.time() - start
print('Number of resulting ligands: ', len(ligand_smiles))
print('Time: ', runtime)
