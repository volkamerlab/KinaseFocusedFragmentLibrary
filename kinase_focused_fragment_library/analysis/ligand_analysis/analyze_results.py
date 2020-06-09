from rdkit import Chem

from .construct_ligand import construct_ligand
from .drug_likeliness import is_drug_like
from .standardize import standardize_mol


def analyze_result(meta, data, original_ligands, chembl):

    try:

        ligand = construct_ligand(meta, data)
        # if ligand could not be constructed, skip
        if not ligand:
            print('ERROR: Ligand could not be constructed:', meta.frag_ids)
            return

        # necessary for Lipinski rule
        # Chem.GetSymmSSSR(ligand)

        ligand = standardize_mol(ligand)
        # if ligand could not be standardized, skip
        if not ligand:
            print('ERROR: Ligand could not be standardized:', meta.frag_ids)
            return

        # Lipinski rule
        lipinski, wt, logp, hbd, hba = is_drug_like(ligand)

        # number of atoms
        n_atoms = ligand.GetNumHeavyAtoms()

        # get INCHI for molecular comparisons
        inchi = Chem.MolToInchi(ligand)

        # search in original ligands
        original = 0
        original_sub = 0
        # exact match in original ligand
        if not original_ligands[original_ligands.inchi == inchi].empty:
            original = 1
        # true substructure of original ligands?
        #elif not original_ligands[original_ligands.mol >= ligand].empty:
        #    original_sub = 1

        # exact chembl match
        chembl_match = 0
        chembl_matches = chembl[chembl.standard_inchi == inchi]
        if not chembl_matches.empty:
            chembl_match = 1

        # save results to dictionary
        ligand_dict = {
            'bonds': [list(i) for i in meta.bonds],
            'frag_ids': list(meta.frag_ids),
            'hba': hba,
            'hbd': hbd,
            'mwt': wt,
            'logp': logp,
            'n_atoms': n_atoms,
            'chembl_match': chembl_match,
            'original': original,
            'original_sub': original_sub
        }

        return ligand_dict

    except Exception as e:

        print(e)
        print('ERROR:', meta.frag_ids)
        return
