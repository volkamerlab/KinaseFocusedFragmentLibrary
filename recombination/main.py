from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from collections import deque  # queue
import time
import sys
from pathlib import Path
import pickle

from metaClasses import Combination, PermutationStep, Fragment, Compound, Port
from get_tuple import get_tuple
from pickle_loader import pickle_loader
from results import results_to_file, add_to_results

start = time.time()

if len(sys.argv) > 1:
    in_arg = int(sys.argv[1])
else:
    in_arg = 5000

path = Path('./tmp')
tmp_files = list(path.glob('tmp_queue*'))
for tmp_file in tmp_files:
    Path.unlink(tmp_file)

path = Path('./results')
out_files = list(path.glob('results*'))
for out_file in out_files:
    Path.unlink(out_file)

# ============================= READ DATA ===============================================

path_to_library = Path('../FragmentLibrary')

# list of folders for each subpocket
folders = list(path_to_library.glob('*'))
subpockets = [str(folder)[-2:] for folder in folders]

# read data

# create dictionary with all fragments for each subpocket
# iterate over all fragments and add each fragmentation site to queue

results = set()  # result set (Combinations)
queue = deque()  # queue containing fragmentation sites to be processed (PermutationSteps containing Compounds)
frags_in_queue = set()  # set containing all fragments that have once been in the queue (Combinations)
frag_set = set()  # only used in initialization for avoiding duplicates in fragment data set (smiles & dummy atoms)

data = {}  # (Fragments)
for folder, subpocket in zip(folders, subpockets):

    file = folder / (subpocket + '.sdf')

    # read molecules
    # keep hydrogen atoms
    suppl = Chem.SDMolSupplier(str(file), removeHs=False)
    mols = [f for f in suppl][:in_arg]

    fragments = []
    for i, fragment in enumerate(mols):

        # ========================== INITIALIZATION ===============================

        fragment = Chem.RemoveHs(fragment)
        frag_id = subpocket + '_' + str(i)

        # store unique atom identifiers
        for a, atom in enumerate(fragment.GetAtoms()):
            frag_atom_id = subpocket + '_' + str(a)
            atom.SetProp('frag_atom_id', frag_atom_id)

        # add all dummy atoms of this fragment to the queue if it has not been in there yet
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if not dummy_atoms:
            continue
        frag_smiles, dummy_set = get_tuple(fragment, dummy_atoms)
        if (frag_smiles, dummy_set) in frag_set:
            continue
        frag_set.add((frag_smiles, dummy_set))

        combo = Combination(frag_ids=frozenset([frag_id]))
        frags_in_queue.add(combo)

        ports = [Port(atom_id=dummy.GetProp('frag_atom_id'), subpocket=subpocket, neighboring_subpocket=dummy.GetProp('subpocket'),
                      bond_type=fragment.GetBondBetweenAtoms(dummy.GetIdx(), dummy.GetNeighbors()[0].GetIdx()).GetBondType())
                 for dummy in dummy_atoms]

        compound = Compound(frag_ids=[frag_id], subpockets=[subpocket], ports=ports, bonds=[])
        for dummy in dummy_atoms:
            ps = PermutationStep(mol=compound, dummy=dummy.GetProp('frag_atom_id'), subpocket=subpocket,
                                 neighboring_subpocket=dummy.GetProp('subpocket'))
            queue.append(ps)

        # store fragment in constant data set
        fragment = Fragment(frag_id=frag_id, subpocket=subpocket, ports=ports)
        fragments.append(fragment)

    data[subpocket] = fragments

n_frags = len(frags_in_queue)

print('Number of fragments: ', n_frags)
print('Number of fragmentation sites: ', len(queue))


# ============================= PERMUTATION ===============================================

# IDEA for frags_in_queue AND results:
# store in file when certain number is reached
# when comparing: load files one by one to compare

count_iterations = 0
count_results = 0
n_tmp_file_out = 0
n_tmp_file_in = 0
n_results_out = 0
limit_q = 100000
limit_r = 100000
n_out = int(limit_q / 2)

results_temp = set()

# while queue not empty
while queue:

    # first element in queue of fragmentation sites to be processed
    # print(len(queue))
    ps = queue.popleft()

    # ========================== TEMP INPUT ================================

    # read back from tmp queue output file if queue is empty
    tmp_q_path = Path('tmp/tmp_queue'+str(n_tmp_file_in)+'.pickle')
    if len(queue) == 0 and tmp_q_path.exists():
        pickle_in = tmp_q_path.open('rb')
        for q_object in pickle_loader(pickle_in):
            queue.append(q_object)
        print('Read ' + str(n_out) + ' queue objects from tmp file', n_tmp_file_in)
        print('Size of queue:', len(queue))
        pickle_in.close()
        Path.unlink(tmp_q_path)
        n_tmp_file_in += 1

    # ========================== TEMP OUTPUT ===============================

    # if queue has reached limit length write part of it to temp output file:
    elif len(queue) >= limit_q:
        tmp_q_path = Path('tmp/tmp_queue'+str(n_tmp_file_out)+'.pickle')
        pickle_out = tmp_q_path.open('wb')
        print('Write ' + str(n_out) + ' queue objects to tmp file', n_tmp_file_out)
        for i in range(n_out):
            ps = queue.pop()  # last element of queue
            pickle.dump(ps, pickle_out)
        print('Size of queue:', len(queue))
        pickle_out.close()
        n_tmp_file_out += 1

    # ======================================================================

    compound = ps.compound
    # dummy atom which is supposed to be replaced with new fragment
    dummy_atom = ps.dummy
    # subpocket to be attached
    neighboring_subpocket = ps.neighboring_subpocket
    # subpocket of current fragment
    subpocket = ps.subpocket

    something_added = False

    # check if subpocket already targeted by this fragment
    if neighboring_subpocket in compound.subpockets:
        # store fragment as ligand if no other open fragmentation sites left
        if len(compound.subpockets) > 1 >= len(compound.ports):
            count_iterations += 1
            combo = Combination(frag_ids=frozenset(compound.frag_ids), bonds=frozenset(compound.bonds))
            # add new result to results
            results_temp.add(combo)
            if len(results_temp) >= limit_r:
                add_to_results(results_temp, results, n_results_out)
                results_temp = set()
            if len(results) >= limit_r:
                count_results += len(results)
                n_results_out = results_to_file(results, n_results_out)
                results = set()
        continue

    # ========================== ITERATION OVER FRAGMENTS ===============================

    # iterate over fragments that might be attached at the current position
    for fragment in data[neighboring_subpocket]:

        # check if fragment has matching fragmentation site
        fragment_ports = [port for port in fragment.ports if port.neighboring_subpocket == subpocket]
        if not fragment_ports:
            continue

        fragment_port = fragment_ports[0]
        compound_port = [port for port in compound.ports if port.atom_id == dummy_atom][0]
        if fragment_port.bond_type != compound_port.bond_type:
            continue

        dummy_atom_2 = fragment_port.atom_id

        # combine fragments
        frag_ids = compound.frag_ids + [fragment.frag_id]
        subpockets = compound.subpockets + [neighboring_subpocket]
        bonds = compound.bonds + [frozenset((dummy_atom, dummy_atom_2))]

        # remove dummy atoms
        ports = [port for port in compound.ports if port.atom_id != dummy_atom]
        ports.extend([port for port in fragment.ports if port.atom_id != dummy_atom_2])

        # ========================== ADD TO RESULTS ===============================

        # if no ports present, ligand is finished
        # if max depth is reached, ligand is also finished, because all subpockets are explored
        combo = Combination(frag_ids=frozenset(frag_ids), bonds=frozenset(bonds))
        if len(ports) == 0 or len(subpockets) == 6:
            count_iterations += 1
            # add new result to results
            results_temp.add(combo)
            if len(results_temp) >= limit_r:
                add_to_results(results_temp, results, n_results_out)
                results_temp = set()
            if len(results) >= limit_r:
                count_results += len(results)
                n_results_out = results_to_file(results, n_results_out)
                results = set()
            something_added = True
            continue

        # ========================== ADD TO QUEUE ===============================

        # check if new molecule is already in queue
        if combo in frags_in_queue:
            continue

        # else add ports of new molecule to queue
        frags_in_queue.add(combo)
        for port in ports:
            new_compound = Compound(frag_ids=frag_ids, subpockets=subpockets, ports=ports, bonds=bonds)
            new_ps = PermutationStep(mol=new_compound, dummy=port.atom_id, subpocket=port.subpocket,
                                     neighboring_subpocket=port.neighboring_subpocket)
            queue.append(new_ps)
        something_added = True

    # ===========================================================================

    # if nothing was added to ps.fragment: store fragment itself as ligand (if it has depth>1 and no other dummy atoms and was not yet in queue)
    if not something_added:

        combo = Combination(frag_ids=frozenset(ps.compound.frag_ids), bonds=frozenset(ps.compound.bonds))
        if combo in frags_in_queue:
            continue

        # ========================== ADD TO RESULTS ===============================

        elif len(ps.compound.subpockets) > 1 >= len(ps.compound.ports):
            count_iterations += 1
            # add new result to results
            results_temp.add(combo)
            if len(results_temp) >= limit_r:
                add_to_results(results_temp, results, n_results_out)
                results_temp = set()
            if len(results) >= limit_r:
                count_results += len(results)
                n_results_out = results_to_file(results, n_results_out)
                results = set()

        # ========================== ADD TO QUEUE ===============================

        # if other dummy atoms are present, remove current dummy (as nothing could be attached there) and add fragment to queue
        elif len(ps.compound.ports) > 1:
            new_ports = [port for port in ps.compound.ports if port.atom_id != ps.dummy]
            new_compound = Compound(frag_ids=ps.compound.frag_ids, subpockets=ps.compound.subpockets, ports=new_ports, bonds=ps.compound.bonds)
            for port in new_compound.ports:
                if port.neighboring_subpocket in new_compound.subpockets:
                    continue
                new_ps = PermutationStep(mol=new_compound, dummy=port.atom_id, subpocket=port.subpocket,
                                         neighboring_subpocket=port.neighboring_subpocket)
                queue.append(new_ps)
                something_added = True
            if something_added:
                frags_in_queue.add(combo)

# ============================= OUTPUT ===============================================

# write remaining results to file
add_to_results(results_temp, results, n_results_out)
results_temp = set()
n_results_out = results_to_file(results, n_results_out)
count_results += len(results)
results = set()

runtime = time.time() - start

# print statistics
print('Number of resulting ligands: ', count_results)
print('Number of ligands including duplicates: ', count_iterations)
print('Overall number of fragments in queue: ', len(frags_in_queue))
print('Time: ', runtime)

stat_path = Path('statistics/statistics_' + str(in_arg) + '.txt')
stat_file = stat_path.open('w')
stat_file.write('Fragments ' + str(n_frags))
stat_file.write('\nLigands ' + str(count_results))
stat_file.write('\nLigands2 ' + str(count_iterations))
stat_file.write('\nQFragments ' + str(len(frags_in_queue)))
stat_file.write('\nTime ' + str(runtime))
stat_file.close()
