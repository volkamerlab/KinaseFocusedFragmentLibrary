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

# iterate over all fragments and store each binding site

count_fragments = 0

queue = deque()

for subpocket in data:

    fragments = [frag for frag in data[subpocket]]

    for fragment in fragments[:10]:

        count_fragments += 1

        dummyAtoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']

        for dummy in dummyAtoms:
            ps = PermutationStep(fragment, dummy, 1)
            queue.append(ps)

print('Number of fragments: ', count_fragments)
print('Number of binding sites: ', len(queue))

# ============================= PERMUTATION ===============================================

# while queue not empty
while queue:

    # first element in queue of binding sites to be processed
    ps = queue.popleft()
    # TO DO: find current subpocket!
    # get subpocket to be attached
    neighboringSubpocket = ps.binding_site.GetProp('subpocket')
    print(neighboringSubpocket)
