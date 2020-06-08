from collections import deque  # queue
import time
import argparse
from pathlib import Path
import pickle

from kinase_focused_fragment_library.recombination.classes_meta import \
    Combination, PermutationStep, Fragment, Compound, Port
from kinase_focused_fragment_library.recombination.get_tuple import get_tuple
from kinase_focused_fragment_library.recombination.pickle_loader import pickle_loader
from kinase_focused_fragment_library.recombination.process_results import \
    results_to_file, add_to_results, process_result
from kinase_focused_fragment_library.recombination.process_queue import add_to_queue
from kinase_focused_fragment_library.recombination.brics_rules import is_brics_bond

from rdkit import Chem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

start = time.time()


def main():

    # ============================= COMMAND LINE ARGUMENTS ===================================

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
    parser.add_argument('-o', '--combinatoriallibrary', type=str, help='output path', required=True)
    parser.add_argument('-n', '--n_frags', type=int, help='Number of input fragments per subpocket', required=False)
    parser.add_argument('-s', '--subpockets', type=str, help='Start from these subpockets only.', required=False, nargs='+')
    parser.add_argument('-d', '--depth', type=int, help='Maximum number of fragments per ligand', default=6, choices=range(2, 7))
    args = parser.parse_args()
    max_depth = args.depth

    output_path = Path(args.combinatoriallibrary)

    # create tmp path
    path = output_path / 'tmp'
    if not path.exists():
        Path.mkdir(path)
    # delete tmp files if present
    tmp_files = path.glob('tmp_queue*')
    for tmp_file in tmp_files:
        Path.unlink(tmp_file)

    # create stat output path
    path = output_path / 'statistics'
    if not path.exists():
        Path.mkdir(path)

    # create output path
    path = output_path / 'results'
    if not path.exists():
        Path.mkdir(path)
    # delete output files if present
    out_files = path.glob('results*')
    for out_file in out_files:
        Path.unlink(out_file)

    # ============================= READ DATA ===============================================

    path_to_library = Path(args.fragmentlibrary)

    # list of folders for each subpocket
    subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2']
    folders = [path_to_library for subpocket in subpockets]

    if args.subpockets:
        start_subpockets = args.subpockets
    else:
        start_subpockets = subpockets

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

        # read molecules, keep hydrogen atoms
        suppl = Chem.SDMolSupplier(str(file), removeHs=False)
        mols = [f for f in suppl][:args.n_frags]

        fragments = []
        for i, fragment in enumerate(mols):

            # ========================== INITIALIZATION ===============================

            fragment = Chem.RemoveHs(fragment)
            frag_id = f'{subpocket}_{i}'

            # store unique atom identifiers
            for a, atom in enumerate(fragment.GetAtoms()):
                frag_atom_id = f'{subpocket}_{a}'
                atom.SetProp('frag_atom_id', frag_atom_id)

            # get all dummy atoms of this fragment except the ones corresponding to the X pool
            dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*' and not a.GetProp('subpocket').startswith('X')]
            if not dummy_atoms:
                continue
            frag_smiles, dummy_set = get_tuple(fragment, dummy_atoms)
            # check if this exact fragment has already been found
            if (frag_smiles, dummy_set) in frag_set:
                continue
            # if not, add this fragment to set of fragments
            frag_set.add((frag_smiles, dummy_set))

            # create dummy atom objects
            ports = [Port(atom_id=dummy.GetProp('frag_atom_id'), subpocket=subpocket, neighboring_subpocket=dummy.GetProp('subpocket'),
                          bond_type=fragment.GetBondBetweenAtoms(dummy.GetIdx(), dummy.GetNeighbors()[0].GetIdx()).GetBondType(),
                          environment=dummy.GetNeighbors()[0].GetProp('environment'))
                     for dummy in dummy_atoms]

            # add all dummy atoms of this fragment to the queue if the fragment lies in a starting pocket
            if subpocket in start_subpockets:

                compound = Compound(frag_ids=[frag_id], subpockets=[subpocket], ports=ports, bonds=[])
                for dummy in dummy_atoms:
                    ps = PermutationStep(mol=compound, dummy=dummy.GetProp('frag_atom_id'), subpocket=subpocket,
                                         neighboring_subpocket=dummy.GetProp('subpocket'))
                    queue.append(ps)

                combo = Combination(frag_ids=frozenset([frag_id]))
                frags_in_queue.add(combo)

            # store fragment in constant data set
            fragment = Fragment(frag_id=frag_id, subpocket=subpocket, ports=ports)
            fragments.append(fragment)

        data[subpocket] = fragments

    n_frags = len(frag_set)

    print('Number of fragments: ', n_frags)
    # print('Number of starting fragmentation sites: ', len(queue))


    # ============================= PERMUTATION ===============================================

    count_iterations = 0
    count_results = 0
    n_tmp_file_out = 0
    n_tmp_file_in = 0
    n_results_out = 0
    limit_q = 500000
    limit_r = 500000
    n_out = int(limit_q / 2)

    results_temp = set()

    # while queue not empty
    while queue:

        # first element in queue of fragmentation sites to be processed
        ps = queue.popleft()

        # ========================== TEMP INPUT ================================

        # read back from tmp queue output file if queue is empty
        tmp_q_path = output_path / ('tmp/tmp_queue'+str(n_tmp_file_in)+'.pickle')
        if len(queue) == 0 and tmp_q_path.exists():
            with open(tmp_q_path, 'rb') as pickle_in:
                for q_object in pickle_loader(pickle_in):
                    queue.append(q_object)
                print('Read ' + str(n_out) + ' queue objects from tmp file', n_tmp_file_in)
                print('Size of queue:', len(queue))

            Path.unlink(tmp_q_path)
            n_tmp_file_in += 1

        # ========================== TEMP OUTPUT ===============================

        # if queue has reached limit length write part of it to temp output file:
        elif len(queue) >= limit_q:
            tmp_q_path = output_path / ('tmp/tmp_queue'+str(n_tmp_file_out)+'.pickle')
            with open(tmp_q_path, 'wb') as pickle_out:
                print('Write ' + str(n_out) + ' queue objects to tmp file', n_tmp_file_out)
                for i in range(n_out):
                    ps = queue.pop()  # last element of queue
                    pickle.dump(ps, pickle_out)
                print('Size of queue:', len(queue))
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

        # ========================== ITERATION OVER FRAGMENTS ===============================

        # iterate over fragments that might be attached at the current position
        for fragment in data[neighboring_subpocket]:

            # check if fragment has matching fragmentation site
            fragment_ports = [port for port in fragment.ports if port.neighboring_subpocket == subpocket]
            if not fragment_ports:
                continue

            fragment_port = fragment_ports[0]
            compound_port = next(port for port in compound.ports if port.atom_id == dummy_atom)
            # check bond types
            if fragment_port.bond_type != compound_port.bond_type:
                continue

            # check environment types
            if not is_brics_bond(fragment_port.environment, compound_port.environment):
                continue

            dummy_atom_2 = fragment_port.atom_id

            # combine fragments
            frag_ids = compound.frag_ids + [fragment.frag_id]
            subpockets = compound.subpockets + [neighboring_subpocket]
            bonds = compound.bonds + [frozenset((dummy_atom, dummy_atom_2))]

            # remove dummy atoms
            ports = [port for port in compound.ports if port.atom_id != dummy_atom]
            ports.extend([port for port in fragment.ports if port.atom_id != dummy_atom_2])

            something_added = True

            # ========================== ADD TO RESULTS ===============================

            # if no ports present, ligand is finished
            # if max depth is reached, ligand is also finished
            combo = Combination(frag_ids=frozenset(frag_ids), bonds=frozenset(bonds))
            if len(ports) == 0 or len(subpockets) == max_depth:
                count_iterations += 1
                # add new result to results
                results, results_temp, n_results_out, count_results = process_result(combo, results_temp, results, limit_r,
                                                                                     n_results_out, count_results, output_path)
                continue

            # ========================== ADD TO QUEUE ===============================

            # check if new molecule is already in queue
            if combo in frags_in_queue:
                continue

            # else add ports of new molecule to queue
            new_compound = Compound(frag_ids=frag_ids, subpockets=subpockets, ports=ports, bonds=bonds)
            added_to_queue = add_to_queue(queue, new_compound)
            if added_to_queue:
                frags_in_queue.add(combo)
            # if nothing was added to queue because no new subpockets: add to results
            else:
                results, results_temp, n_results_out, count_results = process_result(combo, results_temp, results, limit_r,
                                                                                     n_results_out, count_results, output_path)

        # ===========================================================================

        # if nothing was added to ps.fragment: store fragment itself as ligand (if it has depth>1 and no other dummy atoms and was not yet in queue)
        if not something_added:

            combo = Combination(frag_ids=frozenset(ps.compound.frag_ids), bonds=frozenset(ps.compound.bonds))

            # ========================== ADD TO QUEUE ===============================

            # if other dummy atoms are present, remove current dummy (as nothing could be attached there) and add fragment to queue
            if len(ps.compound.ports) > 1:
                new_ports = [port for port in ps.compound.ports if port.atom_id != ps.dummy]
                new_compound = Compound(frag_ids=ps.compound.frag_ids, subpockets=ps.compound.subpockets, ports=new_ports, bonds=ps.compound.bonds)

                something_added = add_to_queue(queue, new_compound)

                if something_added:
                    frags_in_queue.add(combo)

            # ========================== ADD TO RESULTS ===============================

            if (len(ps.compound.subpockets) > 1 >= len(ps.compound.ports)) and not something_added:
                count_iterations += 1
                # add new result to results
                results, results_temp, n_results_out, count_results = process_result(combo, results_temp, results, limit_r,
                                                                                     n_results_out, count_results, output_path)

    # ============================= OUTPUT ===============================================

    # write remaining results to file
    add_to_results(results_temp, results, n_results_out, output_path)
    results_temp = set()
    n_results_out = results_to_file(results, n_results_out, output_path)
    count_results += len(results)
    results = set()

    runtime = time.time() - start

    # print statistics
    print('Number of resulting ligands: ', count_results)
    print('Number of ligands including duplicates: ', count_iterations)
    print('Overall number of fragments in queue: ', len(frags_in_queue))
    print('Time: ', runtime)

    stat_path = output_path / ('statistics/statistics_' + str(args.n_frags) + '.txt')

    with open(stat_path, 'w') as stat_file:
        stat_file.write('Fragments ' + str(n_frags))
        stat_file.write('\nMaxFragments ' + str(max_depth))
        stat_file.write('\nStartSubpockets ' + '_'.join(start_subpockets))
        stat_file.write('\nLigands ' + str(count_results))
        stat_file.write('\nLigands2 ' + str(count_iterations))
        stat_file.write('\nQFragments ' + str(len(frags_in_queue)))
        stat_file.write('\nTime ' + str(runtime))


if __name__ == "__main__":
    main()
