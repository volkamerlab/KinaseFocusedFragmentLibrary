from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)
from pathlib import Path
import time
import matplotlib.pyplot as plt
import sys
import math

sys.path.append('../recombination')
sys.path.append('../recombination/construct_ligands')
from construct_ligand import read_fragment_library
from Thread import LigandAnalysisThread

data = read_fragment_library(Path('../FragmentLibrary'))

# =========================================================================

path_to_results = Path('../recombination/results')
in_paths = list(path_to_results.glob('*.pickle'))

combinatorial_library_folder = Path('../CombinatorialLibrary/')
combinatorial_library_file = combinatorial_library_folder / 'combinatorial_library.pickle'

count_pains = 0
lipinski_ligands, filtered_ligands = 0, 0
wt_ligands = 0
logp_ligands = 0
hbd_ligands = 0
hba_ligands = 0
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

# combinatorial_library_file = combinatorial_library_file.open('wb')

start = time.time()

results = []

# # iterate over ligands
# for in_path in in_paths:
#
#     print(str(in_path))
#     with open(in_path, 'rb') as pickle_in:
#         for meta in pickle_loader(pickle_in):
#
#             result = analyze_result(meta, data, pains)
#             if result is not None:
#                 results.append(result)

n_threads = 1
n_files = len(in_paths)
step_size = math.ceil(n_files/n_threads)

threadList = []

for i in range(0, n_files, step_size):

    if i+step_size > n_files:

        files = in_paths[i:]
        print(files)
        thread = LigandAnalysisThread(files, data)
        threadList.append(thread)

    else:
        files = in_paths[i:i+step_size]
        print(files)
        thread = LigandAnalysisThread(files, data)
        threadList.append(thread)

# start parallel processes
for thread in threadList:
    print('start', thread)
    thread.start()
# join parallel processes:
for thread in threadList:
    print('join', thread)
    thread.join()
# join results
for thread in threadList:
    results.extend(thread.results)

count_ligands = len(results)

# combine results
for result in results:

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

    # pains
    count_pains += result.pains

    # number of atoms
    n_atoms[n] = n_atoms[n] + 1 if n in n_atoms else 1


# combinatorial_library_file.close()
filtered_pains = count_ligands - count_pains

runtime = time.time() - start
print('Number of resulting ligands:', count_ligands)
print('Lipinski rule of 5 fulfilled:', lipinski_ligands, lipinski_ligands/count_ligands)
print('No PAINS found:', filtered_pains, filtered_pains/count_ligands)
print('Molecular weight <= 500:', wt_ligands, wt_ligands/count_ligands)
print('LogP <= 5:', logp_ligands, logp_ligands/count_ligands)
print('HB donors <= 5:', hbd_ligands, hbd_ligands/count_ligands)
print('HB acceptors <= 10:', hba_ligands, hba_ligands/count_ligands)
print('Time: ', runtime)


# plot Lipinski rule
rules = [wt_ligands/count_ligands*100, logp_ligands/count_ligands*100, hbd_ligands/count_ligands*100,
         hba_ligands/count_ligands*100, lipinski_ligands/count_ligands*100, filtered_pains/count_ligands*100]
plt.figure()
ax = plt.bar(range(6), rules)
plt.ylabel('# Ligands [%]')
plt.xticks(range(6), [r'MWT $\leq 500$', r'logP $\leq 5$', r'HBD $\leq 5$', r'HBA $\leq 10$', 'Rule of 5', 'No PAINS'])
plt.yticks()
rects = ax.patches
# calculate percentages
labels = [str(round(n, 2))+'%' for n in rules]
for rect, label in zip(rects, labels):
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height + 0.2, label, ha='center', va='bottom')
plt.savefig(combinatorial_library_folder / 'filters.png')


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
    plt.text(rect1.get_x() + rect1.get_width() / 2, height + 0.2, label1, fontsize=8, ha='center', va='bottom')
    height = rect2.get_height()
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
    plt.text(rect1.get_x() + rect1.get_width() / 2, height + 0.2, label1, fontsize=8, ha='center', va='bottom')
    height = rect2.get_height()
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
