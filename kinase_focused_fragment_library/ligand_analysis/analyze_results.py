from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

from .utils import standardize_mol, construct_ligand


def get_ligand_analysis(meta, data, original_ligands, chembl):

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


def is_drug_like(mol):
    """
    Get Lipinski's rule of five criteria for molecule.

    (If used in loop for multiple molecules, it takes about 1s for 2000 molecules.)

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    tuple of int
        Fulfilled criteria (1) or not (0) for Lipinski's rule of five, and its criteria molecule weight, logP, HBD and
        HBA.
    """

    print(mol)

    mol_wt = 1 if Descriptors.ExactMolWt(mol) <= 500 else 0

    logp = 1 if Descriptors.MolLogP(mol) <= 5 else 0

    hbd = 1 if Lipinski.NumHDonors(mol) <= 5 else 0

    hba = 1 if Lipinski.NumHAcceptors(mol) <= 10 else 0

    lipinski = 1 if mol_wt + logp + hbd + hba >= 3 else 0

    return lipinski, mol_wt, logp, hbd, hba