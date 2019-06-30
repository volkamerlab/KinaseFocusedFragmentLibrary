import parmed as pmd
from rdkit import Chem


phosphate = Chem.MolFromSmiles('O[P](O)(O)O')
phosphate2 = Chem.MolFromSmiles('O[P](=O)(O)O')
phosphatep = Chem.MolFromSmiles('O[P+](O)(O)O')
# ribose = Chem.MolFromSmiles('OC1COC(C1O)CO')
ribose = Chem.MolFromSmiles('OC1COCC1O')


def contains_phosphate(mol):
    return mol.HasSubstructMatch(phosphate) \
    or mol.HasSubstructMatch(phosphate2) \
    or mol.HasSubstructMatch(phosphatep)


def contains_ribose(mol):
    return mol.HasSubstructMatch(ribose)


aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                  "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                  "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


def is_covalent(pdb, pdb_id, chain):

    protein_atom_numbers = []

    # struct = pmd.download_PDB(pdb)
    struct = pmd.load_file('PDB/'+pdb+'.cif')
    
    for res in struct.residues:
        if res.chain == chain:
            if res.name in aa:
                protein_atom_numbers.extend([atom.idx for atom in res.atoms])
            if res.name == pdb_id:
                ligand_atoms = res.atoms
                break

    try:
        ligand_atoms
    except NameError:
        print('ERROR: Ligand not found in PDB file: ', pdb, pdb_id)
        return False

    # atom_numbers = [atom.idx for atom in ligand_atoms]

    for atom in ligand_atoms:
        partners = atom.bond_partners
        for partner in partners:
            # if partner.idx not in atom_numbers:
            if partner.idx in protein_atom_numbers:
                print(atom, partner)
                return True

    return False


# print(is_covalent('4cfn', 'JYM', 'A'))
# print(is_covalent('5fed', '5X4', 'A'))
# print(is_covalent('3w2s', 'W2R', 'A'))
# print(is_covalent('1h08', 'BYP', 'A'))
