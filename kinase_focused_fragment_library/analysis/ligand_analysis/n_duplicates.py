from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.PropertyMol import PropertyMol
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
import pickle
import sys
from pathlib import Path
import pandas as pd

from construct_ligand import construct_ligand, read_fragment_library
from standardize import standardize_mol
sys.path.append('../../recombination')
from pickle_loader import pickle_loader

# load fragments (reduced library)
subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']
data = read_fragment_library(Path('/home/paula/Masterarbeit/FragmentLibrary'), subpockets)

# combinatorial library
combinatorial_library_folder = Path('/home/paula/Masterarbeit/CombinatorialLibrary/')
file_name = combinatorial_library_folder / 'combinatorial_library.pickle'

# create empty dataframe to store inchis (length=number of recombinations)
all_inchis = pd.Series('', index=range(15881010))

# load recombinations 
with open(file_name, 'rb') as pickle_file:

    for i, result in enumerate(pickle_loader(pickle_file)):
        
        # print progress
        if i % 100000 == 0:
            print(i)
        
        # construct molecule
        try: 
            ligand = construct_ligand(result.meta, data)           
        except Exception as e:
            print(e)
            print('ERROR:', i, result.meta.frag_ids)
            continue

        # standardize molecule
        try:    
            std_ligand = standardize_mol(ligand)
        except Exception as e:
            print(e)
            print('ERROR:', i, result.meta.frag_ids)
            continue

        # convert mol to inchi    
        try: 
            inchi = Chem.MolToInchi(std_ligand)
        except Exception as e:
            print(e)
            print('ERROR:', i, result.meta.frag_ids)
            continue
            
        # store inchi in dataframe    
        all_inchis[i] = inchi

# count number of unique InChIs in dataframe
print(all_inchis.nunique())
