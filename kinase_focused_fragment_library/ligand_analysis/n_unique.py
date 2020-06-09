"""
Get indices (and total number) of unique recombined molecules.
To be run on a cluster!

At the moment: Needs to be run from folder:
/kinase_focused_fragment_library/analysis/ligand_analysis/
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from rdkit import Chem

from kinase_focused_fragment_library.ligand_analysis import construct_ligand, read_fragment_library
from kinase_focused_fragment_library.ligand_analysis import standardize_mol

"""
At the moment:
We will need to add the path to the Results class because the pickled object 'combinatorial_library' was
constructed before KinaseFocusedFragmentLibrary became a package and does not know where to search for the class.
Responsible file (that needs to be changed to package-style imports):
/kinase_focused_fragment_library/analysis/ligand_analysis/analyze_results.py

Then the recombination would need to be rerun - and then one could remove the sys statements and use this 
path-independent import:
from kinase_focused_fragment_library.recombination.pickle_loader import pickle_loader
"""

sys.path.append('../../recombination')
from pickle_loader import pickle_loader

# Set default pickle properties for rdkit
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

# set IO folders
data_path = Path('/home/dominique/Documents/Projects/KinaseFocusedFragmentLibraryData')
fragment_library_path = data_path / 'FragmentLibrary'
combinatorial_library_path = data_path / 'CombinatorialLibrary'

# set output paths
combinatorial_library_file = combinatorial_library_path / 'combinatorial_library.pickle'
output_unique_molecule_number_file = combinatorial_library_path / 'unique_molecule_number.txt'
output_unique_molecule_ids_file = combinatorial_library_path / 'unique_molecule_ids.txt'

# load fragments (reduced library)
subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']
fragment_library = read_fragment_library(fragment_library_path, subpockets)

# create dictionary for molecule inchis and IDs
molecule_inchis = {'index': [], 'inchi': []}

# load recombinations 
with open(combinatorial_library_file, 'rb') as pickle_file:

    for molecule_index, result in enumerate(pickle_loader(pickle_file)):
        
        # print progress
        if molecule_index % 100000 == 0:
            print(f'Progress: {molecule_index}')
        
        # construct molecule
        try: 
            molecule = construct_ligand(result.meta, fragment_library)
        except Exception as e:
            print(e)
            print(f'Error: Molecule construction failed: '
                  f'Molecule index {molecule_index} with fragment IDs {result.meta.frag_ids}')
            continue

        # standardize molecule
        try:    
            molecule_standardized = standardize_mol(molecule)
        except Exception as e:
            print(e)
            print(f'Error: Molecule standardization failed: '
                  f'Molecule index {molecule_index} with fragment IDs {result.meta.frag_ids}')
            continue

        # convert mol to inchi    
        try: 
            molecule_inchi = Chem.MolToInchi(molecule_standardized)
        except Exception as e:
            print(e)
            print(f'Error: Mol to inchi conversion failed: '
                  f'Molecule index {molecule_index} with fragment IDs {result.meta.frag_ids}')
            continue
            
        # add inchi with molecule ID
        molecule_inchis['index'].append(molecule_index)
        molecule_inchis['inchi'].append(molecule_inchi)

# create series
molecule_inchis_df = pd.Series(molecule_inchis['inchi'], index=molecule_inchis['index'])

# unique molecule ids (drop inchi duplicates)
unique_molecule_ids = list(molecule_inchis_df.drop_duplicates().index)

# write to file: number of unique recombined molecules
with open(output_unique_molecule_number_file, 'w') as f:
    f.write(f'Number of unique recombined molecules: {len(unique_molecule_ids)}')

# write to file: indices of unique recombined molecules
np.savetxt(output_unique_molecule_ids_file, unique_molecule_ids, fmt='%d')
