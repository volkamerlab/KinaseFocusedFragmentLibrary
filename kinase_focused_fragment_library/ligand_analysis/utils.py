from functools import reduce

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.PropertyMol import PropertyMol

RDLogger.DisableLog('rdApp.*')


def read_fragment_library(path_to_library, subpockets):

    """
    Read fragment library

    Parameters
    ----------
    path_to_library: PosixPath
        path to the folder containing the fragment library

    Returns
    -------
    data: dict(list(RDKit Molecule))
        dictionary with list of fragments for each subpocket

    """

    # list of folders for each subpocket
    folders = [path_to_library for subpocket in subpockets]

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

            fragment = PropertyMol(fragment)

            fragments.append(fragment)

        data[subpocket] = fragments

        n_frags = len(fragments)
        print('Number of fragments in', subpocket,  ':', n_frags)

    return data


def read_original_ligands(frag_dict, path_to_klifs):

    print('Read original ligands.')

    kinases_pdbs = set()

    for subpocket in frag_dict:

        for frag in frag_dict[subpocket]:
            kinases_pdbs.add((frag.GetProp('kinase'), frag.GetProp('_Name')))

    inchis = []
    mols = []
    for kinase, pdb in kinases_pdbs:
        f = path_to_klifs / ('HUMAN/' + kinase + '/' + pdb + '/ligand.mol2')
        ligand = Chem.MolFromMol2File(str(f))

        # standardization
        ligand = standardize_mol(ligand)
        # if ligand could not be standardized, skip
        if not ligand:
            print('Ligand could not be standardized: ', pdb)
            return

        mols.append(ligand)
        inchi = Chem.MolToInchi(ligand)
        inchis.append(inchi)

    print('Number of original ligands :', len(inchis))

    ligands = pd.DataFrame(data=inchis, dtype=str, columns=['inchi'])
    # add molecule column
    ligands['mol'] = mols

    return ligands


def read_chembl_ligands(path_to_chembl):
    """
    Read ChEMBL ligands (ChEMBL compound ID and standardized InChI) from file and add fingerprints.

    Parameters
    ----------
    path_to_chembl : pathlib.Path
        Path to standardized ChEMBL data file.

    Returns
    -------
    pandas.DataFrame
        ChEMBL ligands with each the ChEMBL compound ID, standardized InChI, and fingerprint.
    """

    # read chembl ligands from file
    print('Read', path_to_chembl)
    mols = pd.read_csv(path_to_chembl, header=None, names=['chembl_id', 'standard_inchi'])
    print('Number of ChEMBL molecules:', mols.shape[0])

    # generate fingerprint
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    mols['fingerprint'] = mols.standard_inchi.apply(
        lambda x: rdkit_gen.GetFingerprint(Chem.MolFromInchi(x))
    )

    print(f'ChEMBL data columns: {mols.columns}')

    return mols


def construct_ligand(meta, data):

    """
    Construct a ligand by connecting multiple fragments based on a Combination object

    Parameters
    ----------
    meta: Combination object
        Molecule to be constructed
    data: dict(Mol)
        dictionary containing a list of fragments for each subpocket

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
    combo = reduce(Chem.CombineMols, fragments)

    bonds_matching = True
    ed_combo = Chem.EditableMol(combo)
    replaced_dummies = []
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

        replaced_dummies.extend([dummy_1.GetIdx(), dummy_2.GetIdx()])

    # remove replaced dummy atoms
    replaced_dummies.sort(reverse=True)
    for dummy in replaced_dummies:
        ed_combo.RemoveAtom(dummy)

    ligand = ed_combo.GetMol()

    # replace remaining dummy atoms with hydrogens
    du = Chem.MolFromSmiles('*')
    h = Chem.MolFromSmiles('[H]', sanitize=False)
    ligand = AllChem.ReplaceSubstructs(ligand, du, h, replaceAll=True)[0]
    try:
        ligand = Chem.RemoveHs(ligand)
    except ValueError:
        print(Chem.MolToSmiles(ligand))
        return

    # do not construct this ligand if bond types are not matching
    if not bonds_matching:
        return

    # clear properties
    for prop in ligand.GetPropNames():
        ligand.ClearProp(prop)
    for atom in ligand.GetAtoms():
        atom.ClearProp('frag_atom_id')

    return ligand


def standardize_mol(mol):
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.Normalize(mol)
        mol = rdMolStandardize.Reionize(mol)
        u = rdMolStandardize.Uncharger()
        mol = u.uncharge(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        return mol
    except Exception as e:
        print(e)
        return None


def standardize_inchi(input_inchi):
    try:
        mol = Chem.MolFromInchi(input_inchi)
    except:
        print('ERROR in MolFromInchi:', input_inchi)
        return None
    if not mol:
        print('ERROR in MolFromInchi:', input_inchi)
        return None
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.Normalize(mol)
        mol = rdMolStandardize.Reionize(mol)
        u = rdMolStandardize.Uncharger()
        mol = u.uncharge(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except:
        print('ERROR in standardization:', input_inchi)
        return None
    try:
        inchi = Chem.MolToInchi(mol)
    except:
        print('ERROR in MolToInchi:', input_inchi)
        return None
    return inchi