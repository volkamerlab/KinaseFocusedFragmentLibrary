from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.PropertyMol import PropertyMol
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from collections import deque  # queue
import time
import sys
from pathlib import Path
import pickle
sys.path.append("../fragmentation/")

from add_to_ import add_to_results, add_to_queue, get_tuple
from PermutationStep import PermutationStep
from temp_file import queue_to_tmp

start = time.time()

in_arg = int(sys.argv[1])
limit = in_arg*20  # if queue has reached limit, write fragments in queue to file
# limit = 1000
count_iterations = 0


# queue: 2000 objects ~ 1 MB
def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def pickle_loader(pickle_file):
    try:
        while True:
            yield pickle.load(pickle_file)
    except EOFError:
        pass


# ============================= READ DATA ===============================================

path_to_library = Path('../FragmentLibrary')

# list of folders for each subpocket
folders = list(path_to_library.glob('*'))
subpockets = [str(folder)[-2:] for folder in folders]

# read data

# create dictionary with all fragments for each subpocket
# iterate over all fragments and add each fragmentation site to queue

results = set()  # result set
queue = deque()  # queue containing fragmentation sites to be processed
frags_in_queue = set()  # set containing all fragments that have once been in the queue

data = {}
for folder, subpocket in zip(folders, subpockets):

    file = folder / (subpocket + '.sdf')

    # read molecules
    # keep hydrogen atoms
    suppl = Chem.SDMolSupplier(str(file), removeHs=False)
    mols = [f for f in suppl][:in_arg]

    fragments = []
    for fragment in mols:

        # ========================== INITIALIZATION ===============================

        AllChem.Compute2DCoords(fragment, sampleSeed=1)
        # store unique atom identifiers
        for a, atom in enumerate(fragment.GetAtoms()):
            frag_atom_id = subpocket + '_' + str(a)
            atom.SetProp('frag_atom_id', frag_atom_id)
        # add all dummy atoms of this fragment to the queue if it has not been in there yet
        added = add_to_queue(fragment, frags_in_queue, queue, [subpocket], depth=1)
        if added:
            # store fragment in constant data set
            fragments.append(fragment)

    data[subpocket] = fragments

n_frags = len(frags_in_queue)

print('Number of fragments: ', n_frags)
print('Number of fragmentation sites: ', len(queue))


# # ========================== INITIALIZATION ==============================================
#
# # iterate over all fragments and add each fragmentation site to queue
#
# results = set()  # result set
# queue = deque()  # queue containing fragmentation sites to be processed
# frags_in_queue = set()  # set containing all fragments that have once been in the queue [as tuples: (smiles, subpockets of atoms)]
#
# for subpocket, fragments in data.items():
#
#     for fragment in fragments:
#
#         AllChem.Compute2DCoords(fragment, sampleSeed=1)
#         add_to_queue(fragment, frags_in_queue, queue, [subpocket], depth=1)
#
# n_frags = len(frags_in_queue)
#
# print('Number of fragments: ', n_frags)
# print('Number of fragmentation sites: ', len(queue))


# temporary output file for storing part of the queue

# ============================= PERMUTATION ===============================================

count_exceptions = 0
n_tmp_file_out = 0
n_tmp_file_in = 0

# while queue not empty
while queue:

    # first element in queue of fragmentation sites to be processed
    print(len(queue))
    ps = queue.popleft()

    # read back from tmp queue output file if queue is empty
    tmp_q_path = Path('tmp/tmp_queue'+str(n_tmp_file_in)+'.sdf')
    if len(queue) == 0 and tmp_q_path.exists():
        pickle_in = tmp_q_path.open('rb')
        for q_object in pickle_loader(pickle_in):
            queue.append(q_object)
        print('Read queue objects from file.')
        pickle_in.close()
        Path.unlink(tmp_q_path)
        n_tmp_file_in += 1

    # if queue has reached limit length write part of it to temp output file:
    elif len(queue) >= limit:
        tmp_q_path = Path('tmp/tmp_queue'+str(n_tmp_file_out)+'.sdf')
        n_out = int(limit/2)
        pickle_out = tmp_q_path.open('wb')
        print('Write ' + str(n_out) + ' queue objects to file.')
        for i in range(n_out):
            ps = queue.pop()  # last element of queue
            ps.fragment = PropertyMol(ps.fragment)
            pickle.dump(ps, pickle_out)
        pickle_out.close()
        n_tmp_file_out += 1

    fragment = ps.fragment
    # dummy atom which is supposed to be replaced with new fragment
    dummy_atom = fragment.GetAtomWithIdx(ps.dummy)
    # connecting atom where the new bond will be made
    atom = dummy_atom.GetNeighbors()[0]
    # subpocket to be attached
    neighboring_subpocket = dummy_atom.GetProp('subpocket')
    # subpocket of current fragment
    subpocket = atom.GetProp('subpocket')

    something_added = False

    # check if subpocket already targeted by this fragment
    if neighboring_subpocket in ps.subpockets:
        # store fragment as ligand if no other open fragmentation sites left
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if ps.depth > 1 >= len(dummy_atoms):
            count_iterations += 1
            count = add_to_results(fragment, dummy_atoms, results)
            count_exceptions += count
        continue

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

        dummy_1_id = dummy_atom.GetProp('frag_atom_id')
        dummy_2_id = dummy_atom_2.GetProp('frag_atom_id')

        # combine fragments to one molecule object
        combo = Chem.CombineMols(fragment, fragment_2)

        # find dummy atoms of new combined molecule
        dummy_atom_1 = [a for a in combo.GetAtoms() if a.GetProp('frag_atom_id') == dummy_1_id][0]
        dummy_atom_2 = [a for a in combo.GetAtoms() if a.GetProp('frag_atom_id') == dummy_2_id][0]

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
        # RemoveAtom in editable molecule instead?
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy_atom_1.GetSmarts()))
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy_atom_2.GetSmarts()))

        # skip this fragment if coordinates can not be inferred
        #try:
        #    AllChem.EmbedMolecule(result, randomSeed=1, maxAttempts=1)
        #except Exception:
        #    continue

        # TO DO: Store fragment info
        # Problem: result is stored as smiles in result set -> properties will disappear ...
        # result.SetProp('subpocket', fragment.GetProp('subpocket')+' '+fragment_2.GetProp('subpocket'))
        # result.SetProp('kinase', fragment.GetProp('kinase')+' '+fragment_2.GetProp('kinase'))
        # result.SetProp('complex_pdb', fragment.GetProp('complex_pdb') + ' ' + fragment_2.GetProp('complex_pdb'))
        # print(result.GetProp('subpocket'))

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
        # ps.subpockets.append(neighboring_subpocket)  # ps was changed inplace and thus also in queue !!!
        something_added = add_to_queue(result, frags_in_queue, queue, ps.subpockets+[neighboring_subpocket], depth=ps.depth+1)

    # if nothing was added to ps.fragment: store fragment itself as ligand (if it has depth>1 and no other dummy atoms and was not yet in queue)
    if not something_added:
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if get_tuple(fragment, dummy_atoms) in frags_in_queue:
            continue
        elif ps.depth > 1 >= len(dummy_atoms):
            count_iterations += 1
            count = add_to_results(fragment, dummy_atoms, results)
            count_exceptions += count
        # if other dummy atoms are present, remove current dummy (as nothing could be attached there) and add fragment to queue
        elif len(dummy_atoms) > 1:
            fragment_stripped = Chem.DeleteSubstructs(fragment, Chem.MolFromSmiles(dummy_atom.GetSmarts()))
            add_to_queue(fragment_stripped, frags_in_queue, queue, ps.subpockets, ps.depth)


# ============================= OUTPUT ===============================================

# write statistics to file
runtime = time.time() - start
stat_path = Path('statistics_' + str(in_arg) + '.txt')
stat_file = stat_path.open('w')
stat_file.write('Fragments ' + str(n_frags))
stat_file.write('\nLigands ' + str(len(results)))
stat_file.write('\nLigands2 ' + str(count_iterations))
stat_file.write('\nQFragments ' + str(len(frags_in_queue)))
stat_file.write('\nErrors ' + str(count_exceptions))
stat_file.write('\nTime ' + str(runtime))
stat_file.close()

# print statistics
print('Number of resulting ligands: ', len(results))
print('Number of ligands including duplicates: ', count_iterations)
print('Number of ligands where 3D structure could not be inferred: ', count_exceptions)
print('Overall number of fragments in queue: ', len(frags_in_queue))
print('Time: ', runtime)

# write ligands to file
output_path = Path('../CombinatorialLibrary/combinatorial_library.sdf')
output_file = output_path.open('w')
w = Chem.SDWriter(output_file)
for mol in results:
    w.write(Chem.MolFromSmiles(mol))
w.close()
output_file.close()

# draw some ligands
results = [Chem.RemoveHs(Chem.MolFromSmiles(mol)) for mol in results]
img = Draw.MolsToGridImage(list(results)[:100], molsPerRow=6)
img.save('test.png')

# Problems:

# - Queue and result set explodes from in_arg=10 (52 fragments) on  -> Why??

# - AddBonds does not always actually create a bond -> Output example!

# - Inferring 3D coordinates does not always work -> Output example!

# - Inferring 3D coordinates takes super long! -> Solution: smaller number of attempts

# TO DO:

# Parallelize? Write part of queue temporarily to hard drive?

# remember fragments of which resulting ligand consists