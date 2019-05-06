from rdkit import Chem
from rdkit.Chem import AllChem
# import time


# add finished ligand to results, if ligand can be kekulized
def add_to_results(result, dummy_atoms, results):

    # remove all dummy atoms from finished ligand
    for dummy in dummy_atoms:
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy.GetSmarts()))
    # start = time.time()
    try:
        # infer 3D coordinates
        AllChem.EmbedMolecule(result, randomSeed=1, maxAttempts=1)
    except Exception:
        return 1
    # print('Get 3D coordinates: ', time.time() - start)

    # add molecule to result set as smiles string (to avoid duplicates)
    result = Chem.MolToSmiles(result)
    results.add(result)
    return 0
