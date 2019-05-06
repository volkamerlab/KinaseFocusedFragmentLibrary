from rdkit import Chem
# from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from collections import deque  # queue
import glob
import gc
import time

from PermutationStep import PermutationStep
from add_results import add_to_results

start = time.time()

n_frags = 10

# ============================= READ DATA ===============================================

path_to_library = '../FragmentLibrary'

# list of folders for each subpocket
folders = glob.glob(path_to_library+'/*')
subpockets = [folder[-2:] for folder in folders]

# read data

# create dictionary with SDMolSupplier for each subpocket
data = {}
for i, folder in enumerate(folders):
    subpocket = subpockets[i]
    file = folder + '/' + subpocket + '.sdf'

    # read molecules
    # keep hydrogen atoms
    suppl = Chem.SDMolSupplier(file, removeHs=False)

    data[subpocket] = suppl

# print(data)

# ============================= INITIALIZATION ===============================================

# iterate over all fragments and add each binding site to queue

count_fragments = 0
results = set()  # result set
queue = deque()  # queue containing binding sites to be processed

for subpocket, suppl in data.items():

    for fragment in [f for f in suppl][:n_frags]:
        count_fragments += 1
        # get all binding sites and add to queue
        # if fragment has no open binding sites (fragment = entire ligand), nothing happens
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        for dummy in dummy_atoms:
            ps = PermutationStep(fragment, dummy, 1, subpockets=[subpocket])
            queue.append(ps)

print('Number of fragments: ', count_fragments)
print('Number of binding sites: ', len(queue))

# ============================= PERMUTATION ===============================================

count_exceptions = 0

# while queue not empty
while queue:

    # first element in queue of binding sites to be processed
    ps = queue.popleft()
    print(len(queue))

    fragment = ps.fragment
    # dummy atom which is supposed to be replaced with new fragment
    dummy_atom = ps.binding_site
    # connecting atom where the new bond will be made
    atom = dummy_atom.GetNeighbors()[0]
    # subpocket to be attached
    neighboring_subpocket = dummy_atom.GetProp('subpocket')
    # subpocket of current fragment
    subpocket = atom.GetProp('subpocket')

    # check if subpocket already targeted by this fragment
    if neighboring_subpocket in ps.subpockets:
        # store fragment as ligand if no other open binding sites left
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if ps.depth > 1 >= len(dummy_atoms):  # len(dummy atoms) should always be >= 1 -> change to == 1 ?
            count = add_to_results(fragment, dummy_atoms, results)
            count_exceptions += count
        continue

    something_added = False

    # iterate over fragments that might be attached at the current position
    for fragment_2 in [f for f in data[neighboring_subpocket]][:n_frags]:

        # check if fragment has matching binding site
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

        # find dummy atoms
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
            count = add_to_results(result, dummy_atoms, results)
            count_exceptions += count
            if count == 0:
                something_added = True
            continue

        # else add binding sites of new molecule to queue
        ps.subpockets.append(neighboring_subpocket)
        for dummy in dummy_atoms:
            ps_new = PermutationStep(result, dummy, ps.depth+1, subpockets=ps.subpockets)
            queue.append(ps_new)
            something_added = True

    # if nothing was added to ps.fragment: store fragment itself as ligand (if it has depth>1 and no other dummy atoms)
    if not something_added:
        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        if ps.depth > 1 >= len(dummy_atoms):
            count = add_to_results(fragment, dummy_atoms, results)
            count_exceptions += count


print('Number of resulting ligands: ', len(results))

results = [Chem.RemoveHs(Chem.MolFromSmiles(mol)) for mol in results]
count_unconnected = 0
for mol in results:
    # AllChem.Compute2DCoords(mol)

    smiles = Chem.MolToSmiles(mol)
    if '.' in smiles:
        count_unconnected += 1

img = Draw.MolsToGridImage(list(results)[:100], molsPerRow=6)
img.save('test.png')

print('Number of ligands with unconnected fragments: ', count_unconnected)
print('Number of ligands where 3D structure could not be inferred: ', count_exceptions)

gc.collect()

print('Time: ', time.time() - start)

# Problems:

# - AddBonds does not always actually create a bond

# - Inferring 3D coordinates takes super long! -> Solution: smaller number of attempts
