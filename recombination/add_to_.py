from rdkit import Chem
from rdkit.Chem import AllChem

from PermutationStep import PermutationStep


# add finished ligand to results, if ligand can be kekulized
def add_to_results(result, dummy_atoms, results):

    # TO DO: store fragment PDBs?

    # remove all dummy atoms from finished ligand
    for dummy in dummy_atoms:
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy.GetSmarts()))
    try:
        # infer 3D coordinates
        AllChem.EmbedMolecule(result, randomSeed=1, maxAttempts=1)  # higher number of maxAttempts does not seem to change anything except runtime
    except Exception:
        return 1

    # add molecule to result set as smiles string (to avoid duplicates)
    result = Chem.MolToSmiles(result)
    results.add(result)
    return 0


# add all open fragmentation sites of a fragment to queue, if not already present
# return: boolean: Was fragment added to queue?
def add_to_queue(fragment, frags_in_queue, queue, subpockets, depth):

    dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']

    # if fragment has no open fragmentation sites (fragment = entire ligand), nothing happens
    if not dummy_atoms:
        return False

    # check if fragment is already in queue
    frag_smiles = fragment
    # replace dummys with generic dummys (without atom number)
    for dummy in dummy_atoms:
        frag_smiles = Chem.ReplaceSubstructs(frag_smiles, Chem.MolFromSmiles(dummy.GetSmarts()), Chem.MolFromSmiles('*'))[0]
    frag_smiles = Chem.MolToSmiles(frag_smiles)

    # First approach: almost twice as many ligands as result -> Why??

    #atom_subpockets = [atom.GetProp('subpocket') for atom in fragment.GetAtoms()]

    atom_subpockets = [dummy.GetProp('subpocket') for dummy in dummy_atoms]
    atom_subpockets.append(dummy_atoms[0].GetNeighbors()[0].GetProp('subpocket'))

    # if fragment already in queue, do nothing
    if (frag_smiles, atom_subpockets) in frags_in_queue:
        return False

    # if fragment not yet in queue, add all fragmentation sites of this fragment to queue
    frags_in_queue.append((frag_smiles, atom_subpockets))
    for dummy in dummy_atoms:
        ps_new = PermutationStep(fragment, dummy, depth, subpockets=subpockets)
        queue.append(ps_new)

    return True