from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolDrawing, DrawingOptions

from collections import deque  # queue
import glob

from classes import PermutationStep

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

# iterate over all fragments and store each binding site in queue

count_fragments = 0
results = set()  # result set
queue = deque()  # queue containing binding sites to be processed

for subpocket in data:

    fragments = [frag for frag in data[subpocket]]
    for fragment in fragments[:10]:
        count_fragments += 1

        dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']
        for dummy in dummy_atoms:
            ps = PermutationStep(fragment, dummy, 1, subpockets=[subpocket])
            queue.append(ps)

print('Number of fragments: ', count_fragments)
print('Number of binding sites: ', len(queue))

# ============================= PERMUTATION ===============================================

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

    # iterate over fragments that might be attached at the current position
    for fragment_2 in data[neighboring_subpocket]:

        # check if subpocket already targeted by this fragment
        if neighboring_subpocket in ps.subpockets:
            continue

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

        # remove dummy atoms
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy_atom_1.GetSmarts()))
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy_atom_2.GetSmarts()))
        # fix molecule coordinates (only at the end, after removing dummys from finished ligand?)
        # AllChem.EmbedMolecule(result, randomSeed=1, maxAttempts=0)

        # dummy atoms of new molecule
        dummy_atoms = [a for a in result.GetAtoms() if a.GetSymbol() == '*']

        # if no dummy atoms present, ligand is finished
        if not dummy_atoms:
            # store ligand as result, if molecule can be kekulized
            try:
                # fix coordinates
                AllChem.EmbedMolecule(result, randomSeed=1, maxAttempts=0)
                results.add(result)
            except Exception:
                pass
            continue
        # if max depth is reached, ligand is also finished, because all subpockets are explored
        elif ps.depth + 1 == 6:
            # remove all dummy atoms from finished ligand
            for dummy in dummy_atoms:
                result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy.GetSmarts()))
            try:
                # fix coordinates
                AllChem.EmbedMolecule(result, randomSeed=1, maxAttempts=0)
                results.add(result)
            except Exception:
                pass
            continue

        # else add binding sites of new molecule to queue
        ps.subpockets.append(neighboring_subpocket)
        for dummy in dummy_atoms:
            ps_new = PermutationStep(result, dummy, ps.depth+1, subpockets=ps.subpockets)
            queue.append(ps_new)

print(len(results))
for mol in results:
    AllChem.Compute2DCoords(mol)
img = Draw.MolsToGridImage(list(results)[:20], molsPerRow=2)
img.save('test.png')