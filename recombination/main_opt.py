from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from collections import deque  # queue
import sys
from pathlib import Path
import gc
sys.path.append("../fragmentation/")

from Compounds import Compound, Fragment, Port

in_arg = int(sys.argv[1])

import sys


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


# ======================== READ DATA AND INITIALIZING ===================================

path_to_library = Path('../FragmentLibrary')

# list of folders for each subpocket
folders = list(path_to_library.glob('*'))
subpockets = [str(folder)[-2:] for folder in folders]

queue = deque()  # queue containing fragments to be processed
fragment_set = set()

data = {}
for folder, subpocket in zip(folders, subpockets):

    file = folder / (subpocket + '.sdf')

    # read molecules
    # keep hydrogen atoms
    suppl = Chem.SDMolSupplier(str(file), removeHs=False)
    mols = [f for f in suppl][:in_arg]

    all_fragments = []
    for i, fragment in enumerate(mols):

        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']

        # if fragment has no open fragmentation sites (fragment = entire ligand), nothing happens
        if not dummy_atoms:
            continue

        # check if fragment is already in queue
        frag_smiles = fragment
        # replace dummys with generic dummys (without atom number)
        for dummy in dummy_atoms:
            frag_smiles = Chem.ReplaceSubstructs(frag_smiles, Chem.MolFromSmiles(dummy.GetSmarts()), Chem.MolFromSmiles('*'))[0]
        frag_smiles = Chem.MolToSmiles(frag_smiles)
        # tuple of atom subpockets
        atom_subpockets = tuple([atom.GetProp('subpocket') for atom in fragment.GetAtoms()])

        # if fragment already in queue, do nothing
        if (frag_smiles, atom_subpockets) in fragment_set:
            continue
        # if fragment not yet in queue, add fragment to queue and to data set
        fragment_set.add((frag_smiles, atom_subpockets))

        # create compound object from fragment
        # TO DO: Have unique names in fragment library files! (remove chain etc.?)
        # frag_id = fragment.GetProp('_Name')
        frag_id = i
        ports = []
        for dummy in dummy_atoms:
            neighboring_subpocket = dummy.GetProp('subpocket')
            port = Port(subpocket=subpocket, neighboring_subpocket=neighboring_subpocket)
            # TO DO: How to store atom ids?? They will belong to different fragments!! Do I need them?
            ports.append(port)
        compound = Compound(pdbs=[frag_id], subpockets=[subpocket], ports=ports)

        queue.append(compound)

        # create Fragment object to store in the constant data set with all fragments
        fragment = Fragment(pdb=frag_id, subpocket=subpocket, ports=ports)
        all_fragments.append(fragment)

    data[subpocket] = all_fragments  # constant

# queue now contains fragments instead of fragmentation sites
print('Number of fragments: ', len(queue))

gc.collect()

# ================================= PERMUTATION ========================================

frags_in_queue = set()
results = set()  # result set

# while queue not empty
while queue:

    print(len(queue))

    compound = queue.popleft()

    # queue.reverse()

    for port in compound.ports:

        subpocket = port.subpocket
        neighboring_subpocket = port.neighboring_subpocket

        something_added = False

        # check if neighboring subpocket is already targeted by this fragment
        if neighboring_subpocket in compound.subpockets:
            # add to result if no other port left and compound consists of at least 2 fragments
            if len(compound.pdbs) > 1 and len(compound.ports) == 1:
                result = frozenset([(sp, pdb) for (sp, pdb) in zip(compound.subpockets, compound.pdbs)])
                results.add(result)
            continue

        for fragment in data[neighboring_subpocket]:

            # check if fragment has matching port
            if subpocket not in [p.neighboring_subpocket for p in fragment.ports]:
                continue

            # if compound and fragment have matching ports: Combine
            combo_pdbs = compound.pdbs + [fragment.pdb]
            combo_subpockets = compound.subpockets + [fragment.subpocket]  # fragment.subpocket == neighboring_subpocket
            # remove matching ports
            new_ports = [p for p in fragment.ports if p.neighboring_subpocket != subpocket]
            new_ports.extend([p for p in compound.ports if p.neighboring_subpocket != neighboring_subpocket])
            combo = Compound(pdbs=combo_pdbs, subpockets=combo_subpockets, ports=new_ports)

            # result as set of tuples (set because order does not matter)
            result = frozenset([(sp, pdb) for (sp, pdb) in zip(combo.subpockets, combo.pdbs)])

            # check if result already in queue
            if result in frags_in_queue:
                something_added = True
                continue

            # check if combo should be put to result set
            # if all 6 subpockets are targeted; or if no open ports
            if len(combo.pdbs) > 5 or new_ports == []:
                results.add(result)
                something_added = True
                continue

            # add combo to queue
            frags_in_queue.add(result)
            queue.append(combo)
            something_added = True

        # if nothing was added to compound in this step
        if not something_added:
            result = frozenset([(sp, pdb) for (sp, pdb) in zip(compound.subpockets, compound.pdbs)])
            if result in frags_in_queue:
                continue
            # remove current port from compound, as nothing could be attached there
            new_ports = [p for p in compound.ports if p.neighboring_subpocket != neighboring_subpocket]
            compound.ports = new_ports
            # if compound has at least two fragments and no port is left
            if len(compound.pdbs) > 1 and len(compound.ports) == 0:
                results.add(result)
            # if other ports are present, add new compound without current port to queue
            elif len(compound.ports) > 0:
                new_compound = Compound(pdbs=compound.pdbs, subpockets=compound.subpockets, ports=new_ports)
                queue.append(new_compound)


print('Number of results: ', len(results))
print('Number of combined fragments in queue: ', len(frags_in_queue))


# ================================= CONSTRUCT FRAGMENTS ========================================


data = {}

for folder, subpocket in zip(folders, subpockets):

    file = folder / (subpocket + '.sdf')

    # read molecules
    # keep hydrogen atoms
    suppl = Chem.SDMolSupplier(str(file), removeHs=False)
    mols = [f for f in suppl][:in_arg]

    data[subpocket] = mols

ligands = set()

for result in results:

    print(result)

    all_fragments = []
    for frag in result:
        subpocket = frag[0]
        frag_id = frag[1]
        fragment = data[subpocket][frag_id]
        fragment.SetProp('subpocket', subpocket)
        all_fragments.append(fragment)

    n_frags = len(all_fragments)
    print([f.GetProp('_Name') for f in all_fragments])

    # append first fragment to queue
    fragment = all_fragments[0]
    fragments = [fragment]

    for i in range(n_frags-1):

        fragments_next = []

        for fragment in fragments:

            dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']

            for dummy in dummy_atoms:

                atom = dummy.GetNeighbors()[0]
                subpocket = atom.GetProp('subpocket')  # current subpocket
                neighboring_subpocket = dummy.GetProp('subpocket')
                # fragment in neighboring subpocket (multiple are not possible)
                fragment_2 = [f for f in all_fragments if f.GetProp('subpocket') == neighboring_subpocket]
                # if neighboring subpocket not targeted, move to next dummy
                if not fragment_2:
                    continue
                fragment_2 = fragment_2[0]
                # find matching dummy atom(s)
                dummy_atoms_2 = [a for a in fragment_2.GetAtoms() if a.GetSymbol() == '*' and a.GetProp('subpocket') == subpocket]
                # if there are no matching dummy atoms, move to next dummy
                if not dummy_atoms_2:
                    continue

                # bond type to dummy atom
                bond_type_1 = fragment.GetBondBetweenAtoms(dummy.GetIdx(), atom.GetIdx()).GetBondType()

                for dummy_2 in dummy_atoms_2:

                    # check bond types
                    atom_2 = dummy_2.GetNeighbors()[0]
                    bond_type_2 = fragment_2.GetBondBetweenAtoms(dummy_2.GetIdx(), atom_2.GetIdx()).GetBondType()
                    if bond_type_1 != bond_type_2:
                        continue

                    dummy_smarts = dummy.GetSmarts()
                    dummy_2_smarts = dummy_2.GetSmarts()

                    combo = Chem.CombineMols(fragment, fragment_2)

                    # find dummy atoms of new combined molecule
                    dummy_atom_1 = [a for a in combo.GetAtoms() if a.GetSmarts() == dummy_smarts and a.GetProp('subpocket') == neighboring_subpocket][0]
                    dummy_atom_2 = [a for a in combo.GetAtoms() if a.GetSmarts() == dummy_2_smarts and a.GetProp('subpocket') == subpocket][0]

                    # find atoms to be connected
                    atom_1 = dummy_atom_1.GetNeighbors()[0]
                    atom_2 = dummy_atom_2.GetNeighbors()[0]

                    # add bond between atoms
                    ed_combo = Chem.EditableMol(combo)
                    ed_combo.AddBond(atom_1.GetIdx(), atom_2.GetIdx(), order=bond_type_1)
                    ligand = ed_combo.GetMol()

                    # skip this combo if no bond could be created
                    smiles = Chem.MolToSmiles(ligand)
                    if '.' in smiles:
                        continue

                    # remove dummy atoms
                    ligand = Chem.DeleteSubstructs(ligand, Chem.MolFromSmiles(dummy_atom_1.GetSmarts()))
                    ligand = Chem.DeleteSubstructs(ligand, Chem.MolFromSmiles(dummy_atom_2.GetSmarts()))

                    # # skip this fragment if coordinates can not be inferred
                    # try:
                    #     AllChem.EmbedMolecule(ligand, randomSeed=1, maxAttempts=1)
                    # except Exception:
                    #     continue

                    fragments_next.append(ligand)

        fragments = fragments_next

    print(len(fragments))
    for f in fragments:
        ligands.add(Chem.MolToSmiles(f))

img = Draw.MolsToGridImage([Chem.MolFromSmiles(l) for l in ligands], molsPerRow=3)
img.save('test_opt.png')