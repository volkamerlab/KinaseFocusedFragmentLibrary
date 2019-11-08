import parmed as pmd
from rdkit import Chem


phosphate = Chem.MolFromSmiles('O[P](O)(O)O')
phosphate2 = Chem.MolFromSmiles('O[P](=O)(O)O')
phosphatep = Chem.MolFromSmiles('O[P+](O)(O)O')
ribose = Chem.MolFromSmiles('OC1COCC1O')


def contains_phosphate(mol):
    return mol.HasSubstructMatch(phosphate) \
    or mol.HasSubstructMatch(phosphate2) \
    or mol.HasSubstructMatch(phosphatep)


def contains_ribose(mol):
    return mol.HasSubstructMatch(ribose)


covalent = ['4d9t', '4hct', '4kio']

aa = '\t'.join(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"])


def is_covalent(pdb, pdb_id, chain):

    """
    Check if in a given protein-ligand complex, the ligand is covalently linked to the protein
    by downloading the PDB file and checking the CONECT records

    Parameters
    ----------
    pdb: Str
        PDB ID of a protein-ligand structure
    pdb_id: Str
        PDB ID of the ligand of interest
    chain: Str
        Chain of interest

    Returns
    -------
    True if the ligand is covalently linked to the protein, False otherwise.

    """

    if pdb in covalent:
        return True

    protein_atom_numbers = []

    print('Download PDB', pdb, pdb_id)
    try:
        struct = pmd.download_PDB(pdb)
    except OSError:
        print('ERROR: Could not retrieve PDB', pdb, '\n')
        return False

    # find protein and ligand atoms
    for res in struct.residues:
        if res.chain == chain:
            # protein residue
            if res.name in aa:
                protein_atom_numbers.extend([atom.idx for atom in res.atoms])
            # ligand found
            elif res.name == pdb_id:
                ligand_atoms = res.atoms
                break

    try:
        ligand_atoms
    except NameError:
        print('ERROR: Ligand not found in PDB file: ', pdb, pdb_id, '\n')
        return False

    # check if there is a connection between the ligand and the protein
    for atom in ligand_atoms:
        partners = atom.bond_partners
        for partner in partners:
            if partner.idx in protein_atom_numbers:
                print(atom, partner)
                return True

    return False


# not covalent
# print(is_covalent('3ggf', 'GVD', 'A'))
# print(is_covalent('3bhy', '7CP', 'A'))
# print(is_covalent('6hop', 'GJK', 'A'))
# print(is_covalent('1h07', 'MFQ', 'A'))
# print(is_covalent('1h00', 'FCP', 'A'))

# covalent
# print(is_covalent('4d9t', '0JG', 'A'))
# print(is_covalent('4hct', '18R', 'A'))
# print(is_covalent('4kio', 'G5K', 'A'))
# print(is_covalent('5p9l', '7G9', 'A'))
