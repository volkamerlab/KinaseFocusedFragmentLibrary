from construct_ligand import construct_ligand
from drug_likeliness import is_drug_like
from Result import Result
from rdkit import Chem

from standardize import standardize_mol


def analyze_result(meta, data, original_ligands, chembl):

    ligand, pdbs, fragpdbs = construct_ligand(meta, data)
    # if ligand could not be constructed, skip
    if not ligand:
        print('Ligand could not be constructed: ', meta)
        return

    # necessary for Lipinski rule
    # Chem.GetSymmSSSR(ligand)

    ligand = standardize_mol(ligand)
    # if ligand could not be standardized, skip
    if not ligand:
        print('Ligand could not be standardized: ', meta)
        return

    inchi = Chem.MolToInchi(ligand)

    # Lipinski rule
    lipinski, wt, logp, hbd, hba = is_drug_like(ligand)

    # number of atoms
    n = ligand.GetNumHeavyAtoms()

    # search in original ligands
    original = 0
    original_sub = 0

    # exact match in original ligand
    if not original_ligands[original_ligands.inchi == inchi].empty:
        original = 1

    # true substructure of original ligands?
    elif not original_ligands[original_ligands.mol >= ligand].empty:
        original_sub = 1

    # exact chembl match
    chembl_match = 0
    chembl_matches = chembl[chembl.standard_inchi == inchi]
    if not chembl_matches.empty:
        chembl_match = 1

    # construct Result object
    result = Result(meta, lipinski, wt, logp, hbd, hba, n, original, original_sub, chembl_match)

    return result
