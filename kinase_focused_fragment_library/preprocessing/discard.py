import parmed as pmd
# import pypdb
from rdkit import Chem


phosphate = Chem.MolFromSmiles('O[P](O)(O)O')
phosphate2 = Chem.MolFromSmiles('O[P](=O)(O)O')
phosphatep = Chem.MolFromSmiles('O[P+](O)(O)O')
ribose = Chem.MolFromSmiles('OC1COCC1O')


def get_ligand_from_multi_ligands(ligand):

    """
    Given a molecule object including multiple molecules:
     - return None if one of the molecules includes a phosphate or ribose (substrate)
     - else remove the molecules that consist of <= 14 heavy atoms
     - if more than one molecule is left, return None, else return the remaining molecule

    Parameters
    ----------
    ligand: Mol
        RDKit molecule object which should include multiple molecules

    Returns
    -------
    Mol if a single molecule was extracted, None otherwise

    """

    multi_ligands = Chem.GetMolFrags(ligand, asMols=True)
    # do not use structures including substrates
    phosphate_ligands = [l for l in multi_ligands if contains_phosphate(l) or contains_ribose(l)]
    if phosphate_ligands:
        return
    # get only large ligands
    multi_ligands = [l for l in multi_ligands if l.GetNumHeavyAtoms() > 14]
    # if there is more than one large ligand, discard this structure
    if len(multi_ligands) != 1:
        return
    else:
        return multi_ligands[0]


def contains_phosphate(mol):
    return mol.HasSubstructMatch(phosphate) \
    or mol.HasSubstructMatch(phosphate2) \
    or mol.HasSubstructMatch(phosphatep)


def contains_ribose(mol):
    return mol.HasSubstructMatch(ribose)


covalent = ['4d9t', '4hct', '4kio']
not_covalent = ['2clx', '4cfn']

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

    if pdb in not_covalent:
        return False

    # # download pdb file as string
    # print('Download PDB', pdb, pdb_id)
    # try:
    #     pdb_file = pypdb.get_pdb_file(pdb, filetype='pdb', compression=False)
    # except OSError:
    #     print('ERROR: Could not retrieve PDB', pdb, '\n')
    #     return False
    #
    # lines = pdb_file.split('\n')
    # i = 0
    # line = lines[i]
    # while not line.startswith('LINK') and not line.startswith('ATOM'):
    #     i += 1
    #     line = lines[i]
    # while line.startswith('LINK'):
    #     # check LINK entry for covalent bond
    #     # if (aa in line[17:20] and pdb_id in line[46:50]) or (aa in line[16:20] and pdb_id in line[46:50]):
    #     # print(line[17:20])
    #     if (line[17:20] in aa and line[47:50] == pdb_id and line[26] == chain) or (line[17:20] == pdb_id and line[47:50] == aa and line[51] == chain):
    #         return True
    #     i += 1
    #     line = lines[i]

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
