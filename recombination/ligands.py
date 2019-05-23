from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.PropertyMol import PropertyMol
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from pathlib import Path
import pickle
import sys
sys.path.append("../fragmentation/")


# ============================= FUNCTIONS ===============================================

def pickle_loader(pickle_file):

    """
    Load a pickle file with multiple objects

    Parameters
    ----------
    pickle_file: file object
        input binary pickle file

    Returns
    -------
    Generator object with loaded objects

    """

    try:
        while True:
            yield pickle.load(pickle_file)
    except EOFError:
        pass

# ============================= READ FRAGMENT ===============================================

path_to_library = Path('../FragmentLibrary')

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
        frag_id = subpocket + '_' + str(i)
        fragment.SetProp('frag_id', frag_id)

        # store unique atom identifiers
        for a, atom in enumerate(fragment.GetAtoms()):
            frag_atom_id = subpocket + '_' + str(a)
            atom.SetProp('frag_atom_id', frag_atom_id)
            # atom.SetProp('frag_id', frag_id)

        fragments.append(fragment)

    data[subpocket] = fragments

    n_frags = len(fragments)
    print('Number of fragments in '+subpocket+' :', n_frags)


# ============================= LIGAND CONSTRUCTION ============================================

in_path = Path('meta_library.pickle')
pickle_in = in_path.open('rb')

# iterate over ligands
for meta in pickle_loader(pickle_in):

    frag_ids = meta.frag_ids
    bonds = meta.bonds

    fragments = []
    for frag_id in frag_ids:
        subpocket = frag_id[:2]
        idx = int(frag_id[3:])
        fragment = data[subpocket][idx]
        fragments.append(fragment)

    # combine fragments
    i = 0
    fragment = fragments[i]
    for j in range(len(fragments)-1):

        combo = Chem.CombineMols(fragment, fragments[i+1])
        fragment = combo
        i += 1

    #ed_combo = Chem.EditableMol(combo)
    #for bond in bonds:
    #    dummy_1 =
    #    ed_combo.AddBond(atom_1.GetIdx(), atom_2.GetIdx(), order=bond_type_1)


