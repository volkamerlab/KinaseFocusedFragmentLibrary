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

        ligand = standardize_mol(ligand)
        # if ligand could not be standardized, skip
        if not ligand:
            print('ERROR: Ligand could not be standardized:', meta.frag_ids)
            return

        # Lipinski's rule of five
        lipinski, wt, logp, hbd, hba = is_drug_like(ligand)

        # number of atoms
        n_atoms = ligand.GetNumHeavyAtoms()

        # get InCHI for molecular comparisons
        inchi = Chem.MolToInchi(ligand)

        # ligand has exact match in original ligands?
        original_exact_matches = original_ligands[
            original_ligands.inchi == inchi
        ].index.to_list()

        # ligand has substructure match in original ligands?
        original_substructure_matches = original_ligands[
            original_ligands.mol.apply(lambda x: x.HasSubstructMatch(ligand))
        ].index.to_list()

        # ligand has exact match in ChEMBL?
        chembl_exact_matches = chembl[
            chembl.standard_inchi == inchi
        ].index.to_list()

        # save results to dictionary
        ligand_dict = {
            'bond_ids': [list(i) for i in meta.bonds],
            'fragment_ids': list(meta.frag_ids),
            'hba': hba,
            'hbd': hbd,
            'mwt': wt,
            'logp': logp,
            'n_atoms': n_atoms,
            'chembl_exact': chembl_exact_matches,
            'original_exact': original_exact_matches,
            'original_substructure': original_substructure_matches
        }

        return ligand_dict

    except Exception as e:

        print(e)
        print('ERROR:', meta.frag_ids)
        return
