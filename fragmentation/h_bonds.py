from rdkit import Chem
import numpy as np
import itertools


def infer_h_bonds(ligand, pocket, pocket_mol2, ifp, factory):

    """
    Given a ligand, a kinase binding pocket, and its IFP, find the atoms of the ligands that are interacting with the kinase hinge region via H bonds

    Parameters
    ----------
    ligand: RDKit Mol object
    pocket: RDKit Mol object
    pocket_mol2: Pandas DataFrame
        atom block from mol2 file representing the pocket
    ifp: list(int)
        binary list, interaction fingerprint given by KLIFS
    factory: RDKit feature factory
        built using:
        fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    Returns
    -------
    h_bond_atoms: list(int)
        list of atom ids for the atoms interacting with the hinge region

    """

    h_bond_atoms = set()
    ligand_conf = ligand.GetConformer()
    # hinge region residues
    residues = [46, 47, 48]
    hbd_res = []
    hba_res = []

    # get IFP for residues 46, 47, 48
    for res in residues:
        bit = res*7-7
        bits = ifp[bit:bit+7]
        # H bond donors
        if bits[3] == '1':
            hbd_res.append(res)
        # H bond acceptors
        if bits[4] == '1':
            hba_res.append(res)

    if not hbd_res and not hba_res:
        return h_bond_atoms

    # get ligand acceptors and donors
    feats = factory.GetFeaturesForMol(ligand)
    hbds = list(find_hbds(feats))
    hbas = list(find_hbas(feats))

    # calculate distances and infer atoms hydrogen bonding to the hinge region

    # protein as H bond donor:
    for res in hbd_res:

        # get atoms of current residue
        res_atoms = list(pocket_mol2[pocket_mol2.res_id == res].index)
        # get residue as mol object
        res_mol = mol_from_atom_ids(pocket, res_atoms)
        # initialize ring info
        # res_mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(res_mol)

        # get conformer
        res_mol_conf = res_mol.GetConformer()
        # find h bond donors
        feats = factory.GetFeaturesForMol(res_mol)
        hbd_res_atoms = list(find_hbds(feats))
        # print('Protein residue', res, 'donor atoms:', [(atom_id, res_mol.GetAtomWithIdx(atom_id).GetSymbol()) for atom_id in hbd_res_atoms])
        # distances from all residue hb donors to all ligand hb acceptors
        for p_atom in hbd_res_atoms:
            p_pos = res_mol_conf.GetAtomPosition(p_atom)
            for l_atom in hbas:
                l_pos = ligand_conf.GetAtomPosition(l_atom)
                dist = np.linalg.norm(p_pos - l_pos)
                if dist <= 3.5:
                    # print('donor:', res, p_atom, res_mol.GetAtomWithIdx(p_atom).GetSymbol(),
                    #       'acceptor:', l_atom, ligand.GetAtomWithIdx(l_atom).GetSymbol(), dist)
                    h_bond_atoms.add(l_atom)

    # protein as H bond acceptor:
    for res in hba_res:

        # get atoms of current residue
        res_atoms = list(pocket_mol2[pocket_mol2.res_id == res].index)
        # get residue as mol object
        res_mol = mol_from_atom_ids(pocket, res_atoms)
        # initialize ring info
        # res_mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(res_mol)
        # get conformer
        res_mol_conf = res_mol.GetConformer()
        # find h bond donors
        feats = factory.GetFeaturesForMol(res_mol)
        hba_res_atoms = list(find_hbas(feats))
        # print('Protein residue', res, 'acceptor atoms:', [(atom_id, res_mol.GetAtomWithIdx(atom_id).GetSymbol()) for atom_id in hba_res_atoms])
        # distances from all residue hb acceptors to all ligand hb donors
        for p_atom in hba_res_atoms:
            p_pos = res_mol_conf.GetAtomPosition(p_atom)
            for l_atom in hbds:
                l_pos = ligand_conf.GetAtomPosition(l_atom)
                dist = np.linalg.norm(p_pos - l_pos)
                if dist <= 3.5:
                    # print('acceptor:', res, p_atom, res_mol.GetAtomWithIdx(p_atom).GetSymbol(),
                    #       'donor:', l_atom, ligand.GetAtomWithIdx(l_atom).GetSymbol(), dist)
                    h_bond_atoms.add(l_atom)

    return h_bond_atoms


def mol_from_atom_ids(mol, atom_path):

    """
    Given a molecule and a list of atom numbers, construct the substructure containing these atoms

    Parameters
    ----------
    mol: RDKit Mol object
    atom_path: list(int)
        list of atom IDs

    Returns
    -------
    RDKit Mol object

    """

    bonds = []
    atommap = {}

    for i, j in itertools.combinations(atom_path, 2):
        b = mol.GetBondBetweenAtoms(i, j)
        if b:
            bonds.append(b.GetIdx())

    return Chem.PathToSubmol(mol, bonds, atomMap=atommap)


def find_hbds(feats):

    """
    Find possible hydrogen bond donors

    Parameters
    ----------
    feats: RDKit chemical feature list
        constructed using:
        fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        factory.GetFeaturesForMol(m)

    Returns
    -------
    Generator object containing atom ids for hydrogen bond donor atoms

    """

    for f in feats:
        if f.GetFamily() == 'Donor':
            yield f.GetAtomIds()[0]


def find_hbas(feats):

    """
    Find possible hydrogen acceptors

    Parameters
    ----------
    feats: RDKit chemical feature list
        constructed using:
        fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        factory.GetFeaturesForMol(m)

    Returns
    -------
    Generator object containing atom ids for hydrogen bond acceptor atoms

    """

    for f in feats:
        if f.GetFamily() == 'Acceptor':
            yield f.GetAtomIds()[0]
