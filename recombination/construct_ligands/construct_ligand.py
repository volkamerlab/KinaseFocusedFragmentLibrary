from rdkit import Chem
from pathlib import Path

# ============================= READ FRAGMENT ===============================================

path_to_library = Path('../../FragmentLibrary')

# list of folders for each subpocket
folders = list(path_to_library.glob('*'))
subpockets = [str(folder)[-2:] for folder in folders]

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


# ============================= LIGAND CONSTRUCTION ============================================

def construct_ligand(meta, ligand_smiles=None):

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

    global data

    frag_ids = meta.frag_ids
    bonds = [tuple(bond) for bond in meta.bonds]

    fragments = []
    for frag_id in frag_ids:
        subpocket = frag_id[:2]
        idx = int(frag_id[3:])
        fragment = data[subpocket][idx]
        fragments.append(fragment)

    # combine fragments
    i = 0
    fragment = fragments[i]
    for j in range(len(fragments) - 1):
        combo = Chem.CombineMols(fragment, fragments[i + 1])
        fragment = combo
        i += 1

    # for fragment_a, fragment_b in zip(fragments[:-1], fragments[1:]):
    #     combo = Chem.CombineMols(fragment_a, fragment_b)

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
