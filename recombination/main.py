from rdkit import Chem
from rdkit.Chem import Draw

from collections import deque  # queue
import time
import sys
from pathlib import Path
sys.path.append("../fragmentation/")

from add_to_ import add_to_results, add_to_queue

start = time.time()

input = int(sys.argv[1])
count_iterations = 0

# ============================= READ DATA ===============================================

path_to_library = Path('../FragmentLibrary')

# list of folders for each subpocket
folders = list(path_to_library.glob('*'))
subpockets = [str(folder)[-2:] for folder in folders]

# read data

# create dictionary with all fragments for each subpocket
data = {}
for folder, subpocket in zip(folders, subpockets):

    file = folder / (subpocket + '.sdf')

    # read molecules
    # keep hydrogen atoms
    suppl = Chem.SDMolSupplier(str(file), removeHs=False)
    fragments = [f for f in suppl]

    data[subpocket] = fragments[:input]

# print(data)

# ========================== INITIALIZATION ==============================================

# iterate over all fragments and add each fragmentation site to queue

results = set()  # result set
queue = deque()  # queue containing fragmentation sites to be processed
frags_in_queue = []  # list containing all fragments that have once been in the queue [as tuples: (smiles, subpockets of atoms)]

for subpocket, fragments in data.items():

    for fragment in fragments:

        add_to_queue(fragment, frags_in_queue, queue, [subpocket], depth=0)

n_frags = len(frags_in_queue)

print('Number of fragments: ', n_frags)
print('Number of fragmentation sites: ', len(queue))

stat_file = Path('statistics_' + str(input) + '.txt').open('w')
stat_file.write('Fragments ' + str(n_frags) + '\n')
stat_file.write('Queue ')

# ============================= PERMUTATION ===============================================

count_exceptions = 0

# while queue not empty
while queue:

    # first element in queue of fragmentation sites to be processed
    print(len(queue))
    stat_file.write(str(len(queue)) + ' ')
    ps = queue.popleft()

    fragment = ps.fragment
    # dummy atom which is supposed to be replaced with new fragment
    dummy_atom = ps.dummy
    # connecting atom where the new bond will be made
    atom = dummy_atom.GetNeighbors()[0]
    # subpocket to be attached
    neighboring_subpocket = dummy_atom.GetProp('subpocket')
    # subpocket of current fragment
    subpocket = atom.GetProp('subpocket')

    # check if subpocket already targeted by this fragment
    if neighboring_subpocket in ps.subpockets:
        # store fragment as ligand if no other open fragmentation sites left
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if ps.depth > 1 >= len(dummy_atoms):  # len(dummy atoms) should always be >= 1 -> change to == 1 ?
            count_iterations += 1
            count = add_to_results(fragment, dummy_atoms, results)
            count_exceptions += count
        continue

    something_added = False

    # iterate over fragments that might be attached at the current position
    for fragment_2 in data[neighboring_subpocket]:

        # check if fragment has matching fragmentation site
        dummy_atoms_2 = [a for a in fragment_2.GetAtoms() if a.GetSymbol() == '*'
                         and a.GetProp('subpocket') == subpocket]
        if not dummy_atoms_2:
            continue

        dummy_atom_2 = dummy_atoms_2[0]
        atom_2 = dummy_atom_2.GetNeighbors()[0]

        # check bond types
        bond_type_1 = fragment.GetBondBetweenAtoms(dummy_atom.GetIdx(), atom.GetIdx()).GetBondType()
        bond_type_2 = fragment_2.GetBondBetweenAtoms(dummy_atom_2.GetIdx(), atom_2.GetIdx()).GetBondType()
        if bond_type_1 != bond_type_2:
            continue

        # combine fragments to one molecule object
        combo = Chem.CombineMols(fragment, fragment_2)

        # find dummy atoms of new combined molecule
        dummy_atoms = [a for a in combo.GetAtoms() if a.GetSymbol() == '*']
        dummy_atom_1 = [a for a in dummy_atoms if a.GetProp('subpocket') == neighboring_subpocket][0]
        dummy_atom_2 = [a for a in dummy_atoms if a.GetProp('subpocket') == subpocket][0]
        # find atoms to be connected
        atom_1 = dummy_atom_1.GetNeighbors()[0]
        atom_2 = dummy_atom_2.GetNeighbors()[0]

        # add bond between atoms
        ed_combo = Chem.EditableMol(combo)
        ed_combo.AddBond(atom_1.GetIdx(), atom_2.GetIdx(), order=bond_type_1)
        result = ed_combo.GetMol()

        # skip this fragment if no bond could be created
        smiles = Chem.MolToSmiles(result)
        if '.' in smiles:
            continue

        # remove dummy atoms
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy_atom_1.GetSmarts()))
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy_atom_2.GetSmarts()))

        # dummy atoms of new molecule
        dummy_atoms = [a for a in result.GetAtoms() if a.GetSymbol() == '*']

        # if no dummy atoms present, ligand is finished
        # if max depth is reached, ligand is also finished, because all subpockets are explored
        if (not dummy_atoms) or (ps.depth+1 == len(subpockets)):
            count_iterations += 1
            count = add_to_results(result, dummy_atoms, results)
            count_exceptions += count
            if count == 0:
                something_added = True
            continue

        # else add fragmentation sites of new molecule to queue
        ps.subpockets.append(neighboring_subpocket)
        something_added = add_to_queue(result, frags_in_queue, queue, ps.subpockets, depth=ps.depth+1)

    # if nothing was added to ps.fragment: store fragment itself as ligand (if it has depth>1 and no other dummy atoms)
    if not something_added:
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if ps.depth > 1 >= len(dummy_atoms):
            count_iterations += 1
            count = add_to_results(fragment, dummy_atoms, results)
            count_exceptions += count

# ============================= OUTPUT ===============================================

# write statistics to file
runtime = time.time() - start
stat_file.write('\nLigands ' + str(len(results)))
stat_file.write('\nLigands2 ' + str(count_iterations))
stat_file.write('\nErrors ' + str(count_exceptions))
stat_file.write('\nTime ' + str(runtime))
stat_file.close()

# print statistics
print('Number of resulting ligands: ', len(results))
print('Number of ligands including duplicates: ', count_iterations)
print('Number of ligands where 3D structure could not be inferred: ', count_exceptions)
print('Time: ', runtime)

# write ligands to file
output_path = Path('../CombinatorialLibrary/combinatorial_library.sdf')
output_file = output_path.open('w')
w = Chem.SDWriter(output_file)
for mol in results:
    w.write(Chem.MolFromSmiles(mol))
output_file.close()

# draw some ligands
results = [Chem.RemoveHs(Chem.MolFromSmiles(mol)) for mol in results]

img = Draw.MolsToGridImage(list(results)[:100], molsPerRow=6)
img.save('test.png')

print(len(frags_in_queue))
# Problems:

# - AddBonds does not always actually create a bond -> Output example!

# - Inferring 3D coordinates does not always work -> Output example!

# - Inferring 3D coordinates takes super long! -> Solution: smaller number of attempts

# TO DO:

# Parallelize?

# remember fragments of which resulting ligand consists
