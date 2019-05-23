from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from collections import deque  # queue
import time
import sys
from pathlib import Path
import pickle
sys.path.append("../fragmentation/")

from metaClasses import Combination, PermutationStep, Fragment, Compound, Port
from add_to_ import get_tuple
from pickle_loader import pickle_loader

start = time.time()

if len(sys.argv) > 1:
    in_arg = int(sys.argv[1])
else:
    in_arg = 5000

output_path = Path('meta_library.pickle')
if output_path.exists():
    Path.unlink(output_path)


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
frag_set = set()  # only used in initialization for avoiding duplicates in fragment data set

data = {}
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
            # atom.SetProp('frag_id', frag_id)

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

count_iterations = 0
n_tmp_file_out = 0
n_tmp_file_in = 0
limit = 1000000
n_out = int(limit/2)

# while queue not empty
while queue:

    # first element in queue of fragmentation sites to be processed
    # print(len(queue))
    ps = queue.popleft()

    # ========================== TEMP OUTPUT ===============================

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

    # if queue has reached limit length write part of it to temp output file:
    elif len(queue) >= limit:
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
            results.add(combo)
        continue

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

        # if no ports present, ligand is finished
        # if max depth is reached, ligand is also finished, because all subpockets are explored
        combo = Combination(frag_ids=frozenset(frag_ids), bonds=frozenset(bonds))
        if len(ports) == 0 or len(subpockets) == 6:
            count_iterations += 1
            results.add(combo)
            something_added = True
            continue

        # check if new molecule is already in queue
        new_compound = Compound(frag_ids=frag_ids, subpockets=subpockets, ports=ports, bonds=bonds)
        if new_compound in frags_in_queue:
            continue

        # else add ports of new molecule to queue
        frags_in_queue.add(combo)
        for port in ports:
            new_ps = PermutationStep(mol=new_compound, dummy=port.atom_id, subpocket=port.subpocket,
                                     neighboring_subpocket=port.neighboring_subpocket)
            queue.append(new_ps)
        something_added = True

    # if nothing was added to ps.fragment: store fragment itself as ligand (if it has depth>1 and no other dummy atoms and was not yet in queue)
    if not something_added:

        combo = Combination(frag_ids=frozenset(ps.compound.frag_ids), bonds=frozenset(ps.compound.bonds))
        if combo in frags_in_queue:
            continue

        elif len(ps.compound.subpockets) > 1 >= len(ps.compound.ports):
            count_iterations += 1
            results.add(combo)
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

# ============================= OUTPUT ===============================================

# print statistics
print('Number of resulting ligands: ', len(results))
print('Number of ligands including duplicates: ', count_iterations)
print('Overall number of fragments in queue: ', len(frags_in_queue))

with open(output_path, 'wb') as output_file:
    for result in results:
        pickle.dump(result, output_file)