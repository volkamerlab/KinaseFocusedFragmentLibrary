import parmed as pmd


aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                  "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                  "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


def is_covalent(pdb, pdb_id, chain):

    global aa

    protein_atom_numbers = []

    parm = pmd.download_PDB(pdb)
    for res in parm.residues:
        if res.chain == chain:
            if res.name in aa:
                protein_atom_numbers.extend([atom.idx for atom in res.atoms])
            if res.name == pdb_id:
                ligand_atoms = res.atoms
                break

    try:
        ligand_atoms
    except NameError:
        print('ERROR: Ligand not found in PDB file!\n')
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
