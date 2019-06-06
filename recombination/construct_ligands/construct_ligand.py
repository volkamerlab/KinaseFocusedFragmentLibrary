from rdkit import Chem


# ============================= READ FRAGMENT ===============================================

def read_fragment_library(path_to_library):

    """
    Read fragment library

    Parameters
    ----------
    path_to_library: PosixPath
        path to the folder containing

    Returns
    -------
    data: dict(list(RDKit Molecule))
        dictionary with list of fragments for each subpocket

    """

    # list of folders for each subpocket
    folders = list(path_to_library.glob('*'))
    subpockets = [folder.name for folder in folders]

    data = {}
    for folder, subpocket in zip(folders, subpockets):

        file = folder / (subpocket + '.sdf')

        # read molecules
        # keep hydrogen atoms
        suppl = Chem.SDMolSupplier(str(file), removeHs=False)
        mols = [f for f in suppl]

        fragments = []
        for i, fragment in enumerate(mols):

            fragment = Chem.RemoveHs(fragment)

            # store unique atom identifiers
            for a, atom in enumerate(fragment.GetAtoms()):
                frag_atom_id = f'{subpocket}_{a}'
                atom.SetProp('frag_atom_id', frag_atom_id)
                atom.SetProp('frag_id', fragment.GetProp('complex_pdb'))

            fragments.append(fragment)

        data[subpocket] = fragments

        n_frags = len(fragments)
        print('Number of fragments in', subpocket,  ':', n_frags)

    return data


# ============================= LIGAND CONSTRUCTION ============================================

def construct_ligand(meta, data, ligand_smiles=None):

    """
    Construct a ligand by connecting multiple fragments based on a Combination object

    Parameters
    ----------
    meta: Combination object
        Molecule to be constructed
    ligand_smiles: set(str)
        set of SMILES strings; if the ligand is already in this set, it will not be returned

    Returns
    -------
    ligand: RDKit Molecule or None
        None if the ligand was not constructed
        RDKit Molecule else.

    """

    frag_ids = meta.frag_ids
    bonds = [tuple(bond) for bond in meta.bonds]

    fragments = []
    for frag_id in frag_ids:
        subpocket = frag_id[:2]
        idx = int(frag_id[3:])
        fragment = data[subpocket][idx]
        fragments.append(fragment)

    # combine fragments using map reduce model
    from functools import reduce
    combo = reduce(Chem.CombineMols, fragments)

    bonds_matching = True
    ed_combo = Chem.EditableMol(combo)
    for bond in bonds:

        dummy_1 = next(atom for atom in combo.GetAtoms() if atom.GetProp('frag_atom_id') == bond[0])
        dummy_2 = next(atom for atom in combo.GetAtoms() if atom.GetProp('frag_atom_id') == bond[1])
        atom_1 = dummy_1.GetNeighbors()[0]
        atom_2 = dummy_2.GetNeighbors()[0]

        # check bond types
        bond_type_1 = combo.GetBondBetweenAtoms(dummy_1.GetIdx(), atom_1.GetIdx()).GetBondType()
        bond_type_2 = combo.GetBondBetweenAtoms(dummy_2.GetIdx(), atom_2.GetIdx()).GetBondType()
        if bond_type_1 != bond_type_2:
            bonds_matching = False
            break

        ed_combo.AddBond(atom_1.GetIdx(), atom_2.GetIdx(), order=bond_type_1)

    dummy_atoms = [a.GetIdx() for a in combo.GetAtoms() if a.GetSymbol() == '*']
    dummy_atoms.sort(reverse=True)
    for dummy in dummy_atoms:
        ed_combo.RemoveAtom(dummy)

    ligand = ed_combo.GetMol()

    # do not construct this ligand if bond types are not matching
    if not bonds_matching:
        return

    smiles = Chem.MolToSmiles(ligand)
    if smiles in ligand_smiles:
        return

    # clear properties
    for prop in ligand.GetPropNames():
        ligand.ClearProp(prop)
    for atom in ligand.GetAtoms():
        atom.ClearProp('frag_atom_id')

    return ligand