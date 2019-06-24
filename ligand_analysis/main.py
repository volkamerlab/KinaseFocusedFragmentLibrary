from rdkit import Chem
from rdkit.Chem.PropertyMol import PropertyMol
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)
from pathlib import Path
import time
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import pickle
sys.path.append('../')

from pickle_loader import pickle_loader
from construct_ligand import construct_ligand, read_fragment_library
from drug_likeliness import is_drug_like

data = read_fragment_library(Path('../FragmentLibrary'))

# =========================================================================

path_to_results = Path('../recombination/results')
in_paths = path_to_results.glob('*.pickle')

combinatorial_library_folder = Path('../CombinatorialLibrary/')
combinatorial_library_file = combinatorial_library_folder / 'combinatorial_library.pickle'

count_ligands = 0
lipinski_ligands = 0
wt_ligands = 0
logp_ligands = 0
hbd_ligands = 0
hba_ligands = 0
# with open(out_path, 'wb') as out_file:
# ligand_smiles = set()

subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']

n_per_sp = {}
for subpocket in subpockets:
    n_per_sp[subpocket] = 0

n_sp = {}
for i in range(len(subpockets)):
    n_sp[i+1] = 0

n_atoms = {}

# combinatorial_library_file = combinatorial_library_file.open('wb')

start = time.time()

# iterate over ligands
for in_path in in_paths:

    print(str(in_path))
    with open(in_path, 'rb') as pickle_in:
        for meta in pickle_loader(pickle_in):

            ligand = construct_ligand(meta, data)
            if not ligand:
                continue

            count_ligands += 1

            n_sp[len(meta.frag_ids)] += 1
            for frag_id in meta.frag_ids:
                n_per_sp[frag_id[:2]] += 1

            property_ligand = PropertyMol(ligand)
            # pickle.dump(property_ligand, combinatorial_library_file)
            # ligand_smiles.add(Chem.MolToSmiles(ligand))

            # necessary for Lipinski rule
            # Chem.SanitizeMol(ligand)  # slower
            Chem.GetSymmSSSR(ligand)

            # Lipinski rule
            lipinski, wt, logp, hbd, hba = is_drug_like(ligand)
            if lipinski:
                lipinski_ligands += 1
            wt_ligands += wt
            logp_ligands += logp
            hbd_ligands += hbd
            hba_ligands += hba

            # number of atoms
            n = ligand.GetNumHeavyAtoms()
            if n in n_atoms:
                n_atoms[n] += 1
            else:
                n_atoms[n] = 1

# combinatorial_library_file.close()

runtime = time.time() - start
print('Number of resulting ligands:', count_ligands)
print('Lipinski rule of 5:', lipinski_ligands, lipinski_ligands/count_ligands)
print('Molecular weight <= 500:', wt_ligands, wt_ligands/count_ligands)
print('LogP <= 5:', logp_ligands, logp_ligands/count_ligands)
print('HB donors <= 5:', hbd_ligands, hbd_ligands/count_ligands)
print('HB acceptors <= 10:', hba_ligands, hba_ligands/count_ligands)
print('Time: ', runtime)

# plot Lipinski rule
rules = [wt_ligands/count_ligands*100, logp_ligands/count_ligands*100, hbd_ligands/count_ligands*100,
         hba_ligands/count_ligands*100, lipinski_ligands/count_ligands*100]
plt.figure()
ax = plt.bar(range(5), rules)
plt.ylabel('# Ligands [%]')
plt.xticks(range(5), [r'MWT $\leq 500$', r'logP $\leq 5$', r'HBD $\leq 5$', r'HBA $\leq 10$', 'Lipinski\nrule of 5'])
plt.yticks()
rects = ax.patches
# calculate percentages
labels = [str(round(n/count_ligands*100, 2))+'%' for n in rules]
for rect, label in zip(rects, labels):
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height + 0.2, label, ha='center', va='bottom')
plt.savefig(combinatorial_library_folder / 'lipinski.png')

# plot number of subpockets per ligand
plt.figure()
ax = plt.bar(list(map(int, n_sp.keys())), [n/count_ligands*100 for n in n_sp.values()])
plt.ylabel('# Ligands [%]')
plt.xlabel('# Subpockets')
rects = ax.patches
# calculate percentages
labels = [str(round(n/count_ligands*100, 2))+'%' for n in n_sp.values()]
for rect, label in zip(rects, labels):
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height + 0.2, label, ha='center', va='bottom')
plt.savefig(combinatorial_library_folder / 'n_subpockets.png')

# plot number of fragments/ligands occupying each subpocket
plt.figure()
ax = plt.bar(range(6), [n/count_ligands*100 for n in n_per_sp.values()])
plt.xticks(range(6), n_per_sp.keys())
plt.ylabel('# Ligands [%]')
rects = ax.patches
# calculate percentages
labels = [str(round(n/count_ligands*100, 2))+'%' for n in n_per_sp.values()]
for rect, label in zip(rects, labels):
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height + 0.2, label, ha='center', va='bottom')
plt.savefig(combinatorial_library_folder / 'n_frags.png')

# plot number of atoms per ligand
plt.figure()
plt.bar(n_atoms.keys(), [n/count_ligands*100 for n in n_atoms.values()])
plt.xticks(range(10, max(n_atoms)+1, 5))
plt.xlabel('# Heavy atoms')
plt.ylabel('# Ligands [%]')
plt.savefig(combinatorial_library_folder / 'n_atoms.png')
