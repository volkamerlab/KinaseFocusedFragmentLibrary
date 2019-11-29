from pathlib import Path
import time
import matplotlib.pyplot as plt
import pickle
import pandas as pd

import multiprocessing as mp

import sys
sys.path.append('../../recombination')

from construct_ligand import read_fragment_library
from pickle_loader import pickle_loader
from analyze_results import analyze_result
from novelty import read_original_ligands, read_chembl
import argparse

from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)

# ============================= COMMAND LINE ARGUMENTS ===================================

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
parser.add_argument('-chembl', type=str, help='file with standardized chembl data (InChIs)', required=True)
parser.add_argument('-klifs', type=str, help='path to KLIFS_download folder (original ligands)', required=True)
parser.add_argument('-o', '--combinatoriallibrary', type=str, help='output path', required=True)
args = parser.parse_args()

# ============================= INPUT DATA ================================================

subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']
fragments = read_fragment_library(Path(args.fragmentlibrary), subpockets)

# standardize chembl
# chembl = read_chembl(args.chembl)

# read standardized chembl inchis
print('Read', args.chembl)
chembl = pd.read_csv(args.chembl, header=None, names=['standard_inchi'])
print('Number of ChEMBL molecules:', chembl.shape[0])

# original ligands from KLIFS
path_to_klifs = Path(args.klifs) / 'KLIFS_download'
original_ligands = read_original_ligands(fragments, path_to_klifs)

# output file
combinatorial_library_folder = Path(args.combinatoriallibrary)
combinatorial_library_file = combinatorial_library_folder / 'combinatorial_library.pickle'

# objects create by the recombination algorithm
path_to_results = combinatorial_library_folder / 'results'
in_paths = list(path_to_results.glob('*.pickle'))

# ================================ INITIALIZE =========================================

count_ligands = 0
lipinski_ligands, filtered_ligands = 0, 0
wt_ligands = 0
logp_ligands = 0
hbd_ligands = 0
hba_ligands = 0
originals = 0
original_subs = 0
chembl_match = 0
novel = 0

n_per_sp, n_filtered_per_sp = {}, {}
for subpocket in subpockets:
    n_per_sp[subpocket] = 0
    n_filtered_per_sp[subpocket] = 0

n_sp, n_filtered_sp = {}, {}
for i in range(len(subpockets)):
    n_sp[i+1] = 0
    n_filtered_sp[i+1] = 0

n_atoms = {}
n_atoms_filtered = {}

start = time.time()

metas = []
results = []

n_processes = 2  # mp.cpu_count()
print("Number of processors: ", n_processes)
pool = mp.Pool(n_processes)

# ========================= CONSTRUCT AND ANALYZE LIGANDS ==============================

# iterate over ligands
for in_path in in_paths:

    print(str(in_path))
    with open(str(in_path), 'rb') as pickle_in:

        results.extend(pool.starmap(analyze_result, [(meta, fragments, original_ligands, chembl) for meta in pickle_loader(pickle_in)]))

# ================================ COMBINE RESULTS ======================================

print('Process results.')

combinatorial_library_file = combinatorial_library_file.open('wb')

# combine results
for result in results:

    if result is None:
        continue

    # store in combinatorial library
    pickle.dump(result, combinatorial_library_file)

    count_ligands += 1

    # number of subpockets
    n_sp[result.n_subpockets] += 1
    # occupied subpockets
    for frag_id in result.meta.frag_ids:
        n_per_sp[frag_id[:2]] += 1

    # lipinski rule
    lipinski_ligands += result.lipinski
    wt_ligands += result.mwt
    logp_ligands += result.logp
    hbd_ligands += result.hbd
    hba_ligands += result.hba

    n = result.n_atoms
    # if Lipinski rule fulfilled
    if result.lipinski == 1:
        n_filtered_sp[len(result.meta.frag_ids)] += 1
        n_atoms_filtered[n] = n_atoms_filtered[n] + 1 if n in n_atoms_filtered else 1
        for frag_id in result.meta.frag_ids:
            n_filtered_per_sp[frag_id[:2]] += 1
        filtered_ligands += 1

    # original ligands
    originals += result.original
    original_subs += result.original_sub

    # chembl
    chembl_match += result.chembl_match

    # novel ligand
    if result.chembl_match == 0 and result.original == 0 and result.original_sub == 0:
        novel += 1

    # number of atoms
    n_atoms[n] = n_atoms[n] + 1 if n in n_atoms else 1


# ==================================== OUTPUT ============================================

combinatorial_library_file.close()

runtime = time.time() - start
print('Number of resulting ligands:', count_ligands)
print('Exact match in original ligands:', originals)
print('Substructures of original ligands:', original_subs)
print('Exact match in ChEMBL:', chembl_match)
print('Novel ligand:', novel)
print('Lipinski rule of 5 fulfilled:', lipinski_ligands, lipinski_ligands/count_ligands)
print('Molecular weight <= 500:', wt_ligands, wt_ligands/count_ligands)
print('LogP <= 5:', logp_ligands, logp_ligands/count_ligands)
print('HB donors <= 5:', hbd_ligands, hbd_ligands/count_ligands)
print('HB acceptors <= 10:', hba_ligands, hba_ligands/count_ligands)
print('Time: ', runtime)

# ==================================== PLOTS ============================================

# plot Lipinski rule
rules = [wt_ligands/count_ligands*100, logp_ligands/count_ligands*100, hbd_ligands/count_ligands*100,
         hba_ligands/count_ligands*100, lipinski_ligands/count_ligands*100]
plt.figure()
ax = plt.bar(range(5), rules)
plt.ylabel('# Ligands [%]')
plt.xticks(range(5), [r'MWT $\leq 500$', r'logP $\leq 5$', r'HBD $\leq 5$', r'HBA $\leq 10$', 'Rule of 5'])
plt.yticks()
rects = ax.patches
# calculate percentages
labels = [str(round(n, 2))+'%' for n in rules]
for rect, label in zip(rects, labels):
    height = rect.get_height()
    if height != 100 and height != 0:
        plt.text(rect.get_x() + rect.get_width() / 2, height + 0.2, label, ha='center', va='bottom')
plt.savefig(combinatorial_library_folder / 'lipinski.png')


# plot number of subpockets per ligand
plt.figure()
barWidth, space = 0.45, 0.2
x_values = list(map(int, n_sp.keys()))
ax1 = plt.bar(x_values, [n/count_ligands*100 for n in n_sp.values()], label='All ligands', width=barWidth, edgecolor='white')
x_values = [x+barWidth for x in list(map(int, n_filtered_sp.keys()))]
ax2 = plt.bar(x_values, [n/filtered_ligands*100 for n in n_filtered_sp.values()], label='Filtered ligands',
              width=barWidth, edgecolor='white')
plt.ylabel('# Ligands [%]')
plt.xlabel('# Subpockets')

rects1, rects2 = ax1.patches, ax2.patches
# calculate percentages
labels1 = [str(round(n/count_ligands*100, 1))+'%' for n in n_sp.values()]
labels2 = [str(round(n/filtered_ligands*100, 1))+'%' for n in n_filtered_sp.values()]
for rect1, rect2, label1, label2 in zip(rects1, rects2, labels1, labels2):
    height = rect1.get_height()
    if height != 100 and height != 0:
        plt.text(rect1.get_x() + rect1.get_width() / 2, height + 0.2, label1, fontsize=8, ha='center', va='bottom')
    height = rect2.get_height()
    if height != 100 and height != 0:
        plt.text(rect2.get_x() + rect2.get_width() / 2, height + 0.2, label2, fontsize=8, ha='center', va='bottom')

plt.xticks([r + barWidth/2 for r in range(1, 7)], range(1, 7))
plt.legend(loc='upper left')
plt.savefig(combinatorial_library_folder / 'n_subpockets.png')


# plot number of fragments/ligands occupying each subpocket
plt.figure()
barWidth, space = 0.45, 0.2
ax1 = plt.bar(range(6), [n/count_ligands*100 for n in n_per_sp.values()], label='All ligands', width=barWidth, edgecolor='white')
x_values = [x+barWidth for x in range(6)]
ax2 = plt.bar(x_values, [n/filtered_ligands*100 for n in n_filtered_per_sp.values()], label='Filtered ligands', width=barWidth, edgecolor='white')

plt.ylabel('# Ligands [%]')
rects1, rects2 = ax1.patches, ax2.patches
# calculate percentages
labels1 = [str(round(n/count_ligands*100, 1))+'%' for n in n_per_sp.values()]
labels2 = [str(round(n/filtered_ligands*100, 1))+'%' for n in n_filtered_per_sp.values()]
for rect1, rect2, label1, label2 in zip(rects1, rects2, labels1, labels2):
    height = rect1.get_height()
    if height != 100 and height != 0:
        plt.text(rect1.get_x() + rect1.get_width() / 2, height + 0.2, label1, fontsize=8, ha='center', va='bottom')
    height = rect2.get_height()
    if height != 100 and height != 0:
        plt.text(rect2.get_x() + rect2.get_width() / 2, height + 0.2, label2, fontsize=8, ha='center', va='bottom')

plt.xticks([x + barWidth / 2 for x in range(6)], n_per_sp.keys())
plt.legend(loc='upper right')
plt.savefig(combinatorial_library_folder / 'n_frags.png')


# plot number of atoms per ligand
plt.figure()
x = sorted(n_atoms)
y = [n_atoms[key] for key in sorted(n_atoms)]
plt.plot(x, [n/count_ligands*100 for n in y], label='All ligands')
#plt.bar(n_atoms.keys(), [n/count_ligands*100 for n in n_atoms.values()])
x = sorted(n_atoms_filtered)
y = [n_atoms_filtered[key] for key in sorted(n_atoms_filtered)]
plt.plot(x, [n/filtered_ligands*100 for n in y], label='Filtered ligands')
plt.xticks(range(10, max(n_atoms)+1, 5))
plt.xlabel('# Heavy atoms')
plt.ylabel('# Ligands [%]')
plt.legend(loc='upper right')
plt.savefig(combinatorial_library_folder / 'n_atoms.png')
