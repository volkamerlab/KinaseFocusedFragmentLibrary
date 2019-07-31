from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)
from pathlib import Path
import time
import matplotlib.pyplot as plt
import sys
import pickle

import multiprocessing as mp

sys.path.append('../recombination')
sys.path.append('../recombination/construct_ligands')
from construct_ligand import read_fragment_library
from pickle_loader import pickle_loader
from analyze_results import analyze_result
from novelty import read_inchis, read_scaffolds, read_original_ligands

data = read_fragment_library(Path('../FragmentLibrary'))

original_ligands = read_original_ligands(data)

# chembl = read_chembl('/home/paula/Downloads/chembl_25_chemreps.txt')
chembl = read_inchis('../../data/chembl/chembl.txt')

# kinase inhibitor scaffolds identified by Hu and Bajorath (dx.doi.org/10.1021/jm501237k | J. Med. Chem. 2015, 58, 315−332)
scaffolds = read_scaffolds(['../../data/Kinase_Inhibitors_And_Scaffolds/Ki_Subset/Kinase_Based_Scaffold_Sets_Ki.dat',
                            '../../data/Kinase_Inhibitors_And_Scaffolds/IC50_Subset/Kinase_Based_Scaffold_Sets_IC50.dat'])

# ================================ INITIALIZE =========================================

path_to_results = Path('../recombination/results')
in_paths = list(path_to_results.glob('*.pickle'))

combinatorial_library_folder = Path('../CombinatorialLibrary/')
combinatorial_library_file = combinatorial_library_folder / 'combinatorial_library.pickle'

count_ligands = 0
lipinski_ligands, filtered_ligands = 0, 0
wt_ligands = 0
logp_ligands = 0
hbd_ligands = 0
hba_ligands = 0
originals = 0
original_subs = 0
chembl_match = 0
scaffold = 0
novel = 0
# with open(out_path, 'wb') as out_file:
# ligand_smiles = set()

subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']

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

n_processes = mp.cpu_count()
print("Number of processors: ", n_processes)
pool = mp.Pool(n_processes)

# ========================= CONSTRUCT AND ANALYZE LIGANDS ==============================

# iterate over ligands
for in_path in in_paths:

    print(str(in_path))
    with open(in_path, 'rb') as pickle_in:

        results.extend(pool.starmap(analyze_result, [(meta, data, original_ligands, chembl, scaffolds) for meta in pickle_loader(pickle_in)]))


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

    # scaffold as substructures
    scaffold += result.scaffold

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
# print('Kinase inhibitor scaffold contained:', scaffold)
print('Lipinski rule of 5 fulfilled:', lipinski_ligands, lipinski_ligands/count_ligands)
print('Molecular weight <= 500:', wt_ligands, wt_ligands/count_ligands)
print('LogP <= 5:', logp_ligands, logp_ligands/count_ligands)
print('HB donors <= 5:', hbd_ligands, hbd_ligands/count_ligands)
print('HB acceptors <= 10:', hba_ligands, hba_ligands/count_ligands)
print('Time: ', runtime)


# # plot novelty statistics
# y = [originals/count_ligands*100, original_subs/count_ligands*100, scaffold/count_ligands*100,
#      chembl_match/count_ligands*100, novel/count_ligands*100]
# plt.figure()
# ax = plt.bar(range(5), y)
# plt.ylabel('# Ligands [%]')
# plt.xticks(range(5), ['Original\nligand', 'Substr. of\noriginal ligand', 'Scaffold\ncontained',  'ChEMBL\nmatch', 'Novel\nligand'])
# plt.yticks()
# rects = ax.patches
# # calculate percentages
# labels = [str(round(n, 2))+'%' for n in y]
# for rect, label in zip(rects, labels):
#     height = rect.get_height()
#     if height != 100 and height != 0:
#         plt.text(rect.get_x() + rect.get_width() / 2, height + 0.2, label, ha='center', va='bottom')
# plt.savefig(combinatorial_library_folder / 'novelty.png')


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
