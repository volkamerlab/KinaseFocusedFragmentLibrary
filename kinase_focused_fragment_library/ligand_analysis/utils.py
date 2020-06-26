"""
KinaseFocusedFragmentLibrary
Subpocket-based fragmentation of kinase inhibitors

This module contains utility functions for the analysis of the combinatorial library.
"""

from functools import reduce
import logging

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.PropertyMol import PropertyMol

RDLogger.DisableLog('rdApp.*')
logger = logging.getLogger(__name__)


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
        logger.error(Chem.MolToSmiles(ligand))
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

    logger.info(f'Read fragment library...')
    logger.info(f'Path: {path_to_library}')

    data = {}
    for subpocket in subpockets:

        file = path_to_library / (subpocket + '.sdf')

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
        logger.info(f'Number of fragments in {subpocket}: {n_frags}')

    return data


def read_original_ligands(frag_dict, path_to_klifs, path_combinatorial_library):

    logger.info('Read original ligands...')
    logger.info(f'Path: {path_to_klifs}')

    kinases_pdbs = set()

    for subpocket in frag_dict:

        for frag in frag_dict[subpocket]:
            kinases_pdbs.add((frag.GetProp('kinase'), frag.GetProp('_Name')))

    # save original ligand index to file
    original_ligand_index = pd.DataFrame(kinases_pdbs, columns=['kinase', 'pdb'])
    original_ligand_index.to_csv(path_combinatorial_library / 'original_ligands_index.csv')

    inchis = []
    mols = []
    for kinase, pdb in kinases_pdbs:
        f = path_to_klifs / ('HUMAN/' + kinase + '/' + pdb + '/ligand.mol2')
        ligand = Chem.MolFromMol2File(str(f))

        # standardization
        ligand = standardize_mol(ligand)
        # if ligand could not be standardized, skip
        if not ligand:
            logger.error(f'Ligand could not be standardized: {pdb}')
            return

        mols.append(ligand)
        inchi = Chem.MolToInchi(ligand)
        inchis.append(inchi)

    logger.info(f'Number of original ligands: {len(inchis)}')

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
    logger.info('Read ChEMBL dataset...')
    logger.info(f'Path: {path_to_chembl}')
    mols = pd.read_csv(path_to_chembl)
    logger.info(f'Number of ChEMBL molecules: {mols.shape[0]}')

    # generate fingerprint
    logger.info('Add fingerprints to ChEMBL dataset...')
    mols['fingerprint'] = mols.inchi.apply(get_fingerprint_from_inchi)
    logger.info(f'ChEMBL data columns: {mols.columns.to_list()}')

    # drop rows with any data missing
    mols.dropna(how='any', inplace=True)
    logger.info(f'Number of ChEMBL molecules after dropping empty fingerprints: {mols.shape[0]}')

    return mols


def get_fingerprint_from_inchi(inchi):
    """
    Get fingerprint from InChI.

    Parameters
    ----------
    inchi : str
        InChI.

    Returns
    -------
    rdkit.DataStructs.cDataStructs.ExplicitBitVect
        Fingerprint.
    """

    # get molecule from InChI
    mol = Chem.MolFromInchi(inchi)

    # get fingerprint generator
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)

    # get fingerprint
    if mol is not None:
        return rdkit_gen.GetFingerprint(mol)
    else:
        return None


def standardize_mol(mol):
    """
    Standardize molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    rdkit.Chem.rdchem.Mol or None
        Standardized molecule or None if standardization failed.
    """

    try:

        # sanitize molecule
        Chem.SanitizeMol(mol)

        # remove non-explicit hydrogens
        mol = Chem.RemoveHs(mol)

        # disconnect metals from molecule
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)

        # normalize moleucle
        mol = rdMolStandardize.Normalize(mol)

        # reionize molecule
        mol = rdMolStandardize.Reionize(mol)

        # uncharge molecule (this helps to standardize protonation states)
        u = rdMolStandardize.Uncharger()
        mol = u.uncharge(mol)

        # assign stereochemistry
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

        return mol

    except Exception as e:

        logger.info(f'ERROR in standardization: {e}')
        return None


def convert_mol_to_inchi(mol):
    """
    Convert molecule to InChI.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    str
        InChI.
    """

    try:

        inchi = Chem.MolToInchi(mol)
        return inchi

    except Exception as e:

        logger.info(f'ERROR in MolToInchi: {e}')
        return None


def standardize_smiles_to_inchi(smiles):
    """
    Standardize molecule: input canonical SMILES and output InChI.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.

    Returns
    -------
    str
        Standardized molecule InChI or None if standardization failed.
    """

    # SMILES to ROMol
    molecule = Chem.MolFromSmiles(smiles)

    # standardize ROMol
    molecule_standardized = standardize_mol(molecule)

    # standardized ROMol to InChI
    if molecule_standardized is not None:

        inchi_standardized = convert_mol_to_inchi(molecule_standardized)
        return inchi_standardized

    else:

        return None


def is_drug_like(mol):
    """
    Get Lipinski's rule of five criteria for molecule.

    (If used in loop for multiple molecules, it takes about 1s for 2000 molecules.)

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    tuple of int
        Fulfilled criteria (1) or not (0) for Lipinski's rule of five, and its criteria molecule weight, logP, HBD and
        HBA.
    """

    mol_wt = 1 if Descriptors.ExactMolWt(mol) <= 500 else 0
    logp = 1 if Descriptors.MolLogP(mol) <= 5 else 0
    hbd = 1 if Lipinski.NumHDonors(mol) <= 5 else 0
    hba = 1 if Lipinski.NumHAcceptors(mol) <= 10 else 0
    lipinski = 1 if mol_wt + logp + hbd + hba >= 3 else 0

    return lipinski, mol_wt, logp, hbd, hba
