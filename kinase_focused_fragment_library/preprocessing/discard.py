import parmed as pmd
from rdkit import Chem


PHOSPHATES = [Chem.MolFromSmiles('O[P](O)(O)O'), Chem.MolFromSmiles('O[P](=O)(O)O'), Chem.MolFromSmiles('O[P+](O)(O)O')]
RIBOSE = Chem.MolFromSmiles('OC1COCC1O')

# These ligands are not correctly identified as covalent/non-covalent using the is_covalent function below.
# This problem was identified after correspondence with Albert Kooistra
COVALENT_LIGAND_PDB_IDS = ['4d9t', '4hct', '4kio']
NON_COVALENT_LIGAND_PDB_IDS = ['2clx', '4cfn']

# three letter codes of all amino acids (used to identify protein atoms in PDB file)
AMINO_ACIDS = '\t'.join(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"])


def get_ligand_from_multi_ligands(ligand):

    """
    Given a molecule object including multiple molecules:
     - return None if one of the molecules includes a phosphate or ribose (substrate)
     - else remove the molecules that consist of <= 14 heavy atoms
     - if more than one molecule is left, return None, else return the remaining molecule

    Parameters
    ----------
    ligand: rdkit.Chem.Mol
        RDKit molecule object which should include multiple molecules

    Returns
    -------
    rdkit.Chem.Mol or None
        Molecule if a single molecule was extracted, None otherwise

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
    """
    Check if given molecule contains a phosphate group
    """
    return mol.HasSubstructMatch(PHOSPHATES[0]) \
    or mol.HasSubstructMatch(PHOSPHATES[1]) \
    or mol.HasSubstructMatch(PHOSPHATES[2])


def contains_ribose(mol):
    """
    Check if given molecule contains a ribose
    """
    return mol.HasSubstructMatch(RIBOSE)


def is_covalent(pdb, ligand_pdb, chain):

    """
    Check if in a given protein-ligand complex, the ligand is covalently linked to the protein
    by downloading the PDB file and checking the CONECT records

    Parameters
    ----------
    pdb: str
        PDB ID of a protein-ligand structure
    ligand_pdb: str
        PDB ID of the ligand of interest
    chain: str
        Chain of interest

    Returns
    -------
    True if the ligand is covalently linked to the protein, False otherwise.

    """

    # These ligands are not correctly identified as covalent/non-covalent using this function.
    # This problem was identified after correspondence with Albert Kooistra
    if pdb in COVALENT_LIGAND_PDB_IDS:
        return True
    if pdb in NON_COVALENT_LIGAND_PDB_IDS:
        return False

    protein_atom_numbers = []

    print('Download PDB', pdb, ligand_pdb)
    try:
        struct = pmd.download_PDB(pdb)
    except OSError:
        print('ERROR: Could not retrieve PDB', pdb, '\n')
        return False

    ligand_atoms = None

    # find protein and ligand atoms
    for res in struct.residues:
        if res.chain == chain:
            # protein residue
            if res.name in AMINO_ACIDS:
                protein_atom_numbers.extend([atom.idx for atom in res.atoms])
            # ligand found
            elif res.name == ligand_pdb:
                ligand_atoms = res.atoms
                break

    if ligand_atoms is None:
        print('ERROR: Ligand not found in PDB file: ', pdb, ligand_pdb, '\n')
        return False

    # check if there is a connection between the ligand and the protein
    for atom in ligand_atoms:
        partners = atom.bond_partners
        for partner in partners:
            if partner.idx in protein_atom_numbers:
                print(atom, partner)
                return True

    return False
