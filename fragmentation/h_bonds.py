from rdkit import Chem
import numpy as np
import itertools


def infer_h_bonds(ligand, pocket, pocket_mol2, ifp, factory):

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

    print(hba_res, hbd_res)

    # get ligand acceptors and donors
    feats = factory.GetFeaturesForMol(ligand)
    hbds = list(find_hbds(feats))
    hbas = list(find_hbas(feats))
    print(hbds, hbas)

    # calculate distances and infer atoms hydrogen bonding to the hinge region

    # protein as H bond donor:
    for res in hbd_res:

        # get atoms of current residue
        hbd_res_atoms = list(pocket_mol2[pocket_mol2.res_id == res].atom_id)
        # get residue as mol object
        res_mol = mol_from_atom_ids(pocket, hbd_res_atoms)
        # initialize ring info
        res_mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(res_mol)
        # get conformer
        res_mol_conf = res_mol.GetConformer()
        # find h bond donors
        feats = factory.GetFeaturesForMol(res_mol)
        hbd_res_atoms = find_hbds(feats)
        # distances from all residue hb donors to all ligand hb acceptors
        for p_atom in hbd_res_atoms:
            p_pos = res_mol_conf.GetAtomPosition(p_atom)
            for l_atom in hbas:
                l_pos = ligand_conf.GetAtomPosition(l_atom)
                dist = np.linalg.norm(p_pos - l_pos)
                if dist <= 4.5:
                    print('donor:', res, p_atom, res_mol.GetAtomWithIdx(p_atom).GetSymbol(), 'acceptor:', l_atom, ligand.GetAtomWithIdx(l_atom).GetSymbol(), dist)
                    h_bond_atoms.add(l_atom)

    # protein as H bond acceptor:
    for res in hba_res:

        # get atoms of current residue
        hba_res_atoms = list(pocket_mol2[pocket_mol2.res_id == res].atom_id)
        # get residue as mol object
        res_mol = mol_from_atom_ids(pocket, hba_res_atoms)
        # initialize ring info
        res_mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(res_mol)
        # get conformer
        res_mol_conf = res_mol.GetConformer()
        # find h bond donors
        feats = factory.GetFeaturesForMol(res_mol)
        hba_res_atoms = find_hbas(feats)
        # distances from all residue hb acceptors to all ligand hb donors
        for p_atom in hba_res_atoms:
            p_pos = res_mol_conf.GetAtomPosition(p_atom)
            for l_atom in hbds:
                l_pos = ligand_conf.GetAtomPosition(l_atom)
                dist = np.linalg.norm(p_pos - l_pos)
                if dist <= 4.5:
                    print('acceptor:', res, p_atom, res_mol.GetAtomWithIdx(p_atom).GetSymbol(), 'donor:', l_atom, ligand.GetAtomWithIdx(l_atom).GetSymbol(), dist)
                    h_bond_atoms.add(l_atom)

    return h_bond_atoms


def mol_from_atom_ids(mol, atom_path):

    bonds = []
    atommap = {}

    for i, j in itertools.combinations(atom_path, 2):
        b = mol.GetBondBetweenAtoms(i, j)
        if b:
            bonds.append(b.GetIdx())

    return Chem.PathToSubmol(mol, bonds, atomMap=atommap)


def find_hbds(feats):

    for f in feats:
        if f.GetFamily() == 'Donor':
            yield f.GetAtomIds()[0]


def find_hbas(feats):

    for f in feats:
        if f.GetFamily() == 'Acceptor':
            yield f.GetAtomIds()[0]
