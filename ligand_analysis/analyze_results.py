from construct_ligand import construct_ligand
from drug_likeliness import is_drug_like
from Result import Result

# global pains
from rdkit.Chem.FilterCatalog import *
params = FilterCatalogParams()
# Build a catalog from all PAINS (A, B and C)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains = FilterCatalog(params)


def analyze_result(meta, data, original_ligands):

    global pains

    ligand = construct_ligand(meta, data)
    # if ligand could not be constructed, skip
    if not ligand:
        print('Ligand could not be constructed: ', meta)
        return

    # necessary for Lipinski rule
    Chem.GetSymmSSSR(ligand)

    # Lipinski rule
    lipinski, wt, logp, hbd, hba = is_drug_like(ligand)

    # PAINS substructure search
    match = pains.GetFirstMatch(ligand)

    pains_found = 0 if match is None else 1

    # number of atoms
    n = ligand.GetNumHeavyAtoms()

    smiles = Chem.MolToSmiles(ligand)
    # search in original ligands
    original = 0
    original_sub = 0
    # for l in original_ligands:
    #     s = Chem.MolToSmiles(l)
    #     found exact match
    #     if s == smiles:
    #         original = 1
    #         break
    #     elif ligand.HasSubstructMatch(l):
    #         original_sub = 1
    #         break

    # chembl
    chembl_match = 0

    # construct Result object
    result = Result(meta, lipinski, wt, logp, hbd, hba, pains_found, n, original, original_sub)

    return result
