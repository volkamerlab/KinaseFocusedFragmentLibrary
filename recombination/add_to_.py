from rdkit import Chem
from rdkit.Chem import AllChem

from PermutationStep import PermutationStep


def add_to_results(result, dummy_atoms, results):

    """
    Adds a fragment to the result set after
    - removing all dummy atoms from the fragment
    - inferring coordinates of the fragment

    Parameters
    ----------
    result: RDKit Mol object
        fragment to be added to results
    dummy_atoms: list(Atom)
        dummy Atom objects of this fragment
    results: set(string)
        result set containing SMILES strings

    Returns
    -------
    int
        0 if fragment was added to the result set
        1 if not, because the molecule could not be kekulized (or other exceptions)

    """

    # remove all dummy atoms from finished ligand
    for dummy in dummy_atoms:
        result = Chem.DeleteSubstructs(result, Chem.MolFromSmiles(dummy.GetSmarts()))
    try:
        # infer coordinates
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

    """
    Add all fragmentation sites (dummy atoms) of a fragment to the queue if
    - fragment has dummy atoms
    - fragment is not already in queue (comparison by smiles, dummy atoms and subpockets of atoms)

    Parameters
    ----------
    fragment: RDKit Mol object
        fragment to be added to queue
    frags_in_queue: set(tuple)
        set of tuples:
        first element: SMILES string of a fragment
        second element: frozenset of dummy atom tuples: (frag_atom_id, subpocket)
    queue: deque(PermutationStep)
        queue containing PermutationStep objects
    subpockets:
        list of the subpockets the fragment is targeting
    depth:
        number of original fragments that the fragment consists of

    Returns
    -------
    Boolean
        True if fragment was added to the queue OR is already in the queue
        False if not, because it has no dummy atoms
    Boolean
        True if fragment was added to the queue
        False if not, because it was already in there

    """

    dummy_atoms = [a for a in fragment.GetAtoms() if a.GetSymbol() == '*']

    # if fragment has no open fragmentation sites (fragment = entire ligand), nothing happens
    if not dummy_atoms:
        return False, False

    # check if fragment is already in queue
    frag_smiles, dummy_set = get_tuple(fragment, dummy_atoms)

    # if fragment already in queue, do nothing
    if (frag_smiles, dummy_set) in frags_in_queue:
        return True, False

    # if fragment not yet in queue, add all fragmentation sites of this fragment to queue
    frags_in_queue.add((frag_smiles, dummy_set))
    for dummy in dummy_atoms:
        ps_new = PermutationStep(fragment, dummy.GetIdx(), depth, subpockets=subpockets)
        queue.append(ps_new)

    return True, True


def get_tuple(fragment, dummy_atoms):

    """
    For a given fragment, returns:
    - smiles string with generic dummy atoms (dummy labels removed)
    - dummy atoms as tuples of frag_atom_id and subpocket (of the dummy = neighboring subpocket of the fragment)

    Parameters
    ----------
    fragment: RDKit Mol object
    dummy_atoms: list(RDKit Atom objects)
        list of all dummy atoms of the fragment

    Returns
    -------
    String
        SMILES string of the fragment

    frozenset(tuple)
        frozenset of tuples for each dummy atom containing the frag_atom_id and the subpocket of the dummy
    """

    frag_smiles = fragment
    # replace dummys with generic dummys (without atom number)
    # dummy tuple: (frag_atom_id, neighboring_subpocket), e.g. (AP_4, FP)
    dummy_set = []
    for dummy in dummy_atoms:
        frag_smiles = Chem.ReplaceSubstructs(frag_smiles, Chem.MolFromSmiles(dummy.GetSmarts()), Chem.MolFromSmiles('*'))[0]
        dummy_tuple = (dummy.GetProp('frag_atom_id'), dummy.GetProp('subpocket'))
        dummy_set.append(dummy_tuple)
    frag_smiles = Chem.MolToSmiles(frag_smiles)

    dummy_set = frozenset(dummy_set)

    return frag_smiles, dummy_set
