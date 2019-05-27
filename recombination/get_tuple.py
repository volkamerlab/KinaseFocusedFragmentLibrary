from rdkit import Chem


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