import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def read_fragment_library(path_to_lib):
    """
    Read fragment library from sdf files (one file per subpocket).

    Parameters
    ----------
    path_to_lib : str
        Path to fragment library folder.


    Returns
    -------
    dict of pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, (values) for each subpocket (key).
    """
    # list of folders for each subpocket
    subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2', 'X']

    data = {}

    # iterate over subpockets
    for subpocket in subpockets:

        fragments = _read_subpocket_fragments(subpocket, path_to_lib)

        if fragments is not None:
            data[subpocket] = fragments

        else:
            pass

    return data


def _read_subpocket_fragments(subpocket, path_to_lib):
    """
    Read fragments for input subpocket.

    Parameters
    ----------
    subpocket : str
        Subpocket name, i.e. AP, SE, FP, GA, B1, or B2.
    path_to_lib : str
        Path to fragment library folder.

    Returns
    -------
    pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, for input subpocket.
    """

    try:
        mol_supplier = Chem.SDMolSupplier(str(path_to_lib / f'{subpocket}.sdf'), removeHs=False)
    except OSError:
        return None

    data = []

    for mol in mol_supplier:
        # Replace dummy atoms with hydrogens in fragments
        dummy = Chem.MolFromSmiles('*')
        hydrogen = Chem.MolFromSmiles('[H]', sanitize=False)
        mol_wo_dummy = AllChem.ReplaceSubstructs(mol, dummy, hydrogen, replaceAll=True)[0]

        # Remove all hydrogens but explicit hydrogens
        mol_wo_dummy = Chem.RemoveHs(mol_wo_dummy)
        mol_w_dummy = Chem.RemoveHs(mol)

        # Generate SMILES
        smiles_wo_dummy = Chem.MolToSmiles(mol_wo_dummy)
        smiles_w_dummy = Chem.MolToSmiles(mol_w_dummy)

        # 2D coordinates
        AllChem.Compute2DCoords(mol_wo_dummy)
        AllChem.Compute2DCoords(mol_w_dummy)

        # Add property information stored for each fragment, e.g. kinase group
        data.append(
            [
                mol_wo_dummy,
                mol_w_dummy,
                mol,
                mol_w_dummy.GetProp('kinase'),
                mol_w_dummy.GetProp('family'),
                mol_w_dummy.GetProp('group'),
                mol_w_dummy.GetProp('complex_pdb'),
                mol_w_dummy.GetProp('ligand_pdb'),
                mol_w_dummy.GetProp('alt'),
                mol_w_dummy.GetProp('chain'),
                mol_w_dummy.GetProp('atom.prop.subpocket'),
                mol_w_dummy.GetProp('atom.prop.environment'),
                smiles_wo_dummy,
                smiles_w_dummy
            ]
        )

    fragment_library = pd.DataFrame(
        data,
        columns='ROMol ROMol_dummy ROMol_original kinase family group complex_pdb ligand_pdb alt chain atom_subpockets atom_environments smiles smiles_dummy'.split()
    )
    fragment_library['subpocket'] = subpocket

    return fragment_library