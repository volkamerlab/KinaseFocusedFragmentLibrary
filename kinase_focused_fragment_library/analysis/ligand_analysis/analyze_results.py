from construct_ligand import construct_ligand
from drug_likeliness import is_drug_like
from Result import Result

from standardize import standardize_mol, standardize_smiles

# global pains
from rdkit.Chem.FilterCatalog import *
params = FilterCatalogParams()
# Build a catalog from all PAINS (A, B and C)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains = FilterCatalog(params)


def analyze_result(meta, data, original_ligands, chembl, scaffolds):

    global pains

    ligand, pdbs, fragpdbs = construct_ligand(meta, data)
    # if ligand could not be constructed, skip
    if not ligand:
        print('Ligand could not be constructed: ', meta)
        return

    # TODO: standardize molecules instead of GetSymmSSSR and removing p+1 from InChI

    # necessary for Lipinski rule
    # Chem.GetSymmSSSR(ligand)
    # ligand = standardize_mol(ligand)

    # standardize molecule
    smiles = Chem.MolToSmiles(ligand)
    smiles = standardize_smiles(smiles)
    ligand = Chem.MolFromSmiles(smiles)
    # inchi = Chem.MolToInchi(ligand)  #.replace('/p+1', '')

    # Lipinski rule
    lipinski, wt, logp, hbd, hba = is_drug_like(ligand)

    # PAINS substructure search
    match = pains.GetFirstMatch(ligand)

    pains_found = 0 if match is None else 1

    # number of atoms
    n = ligand.GetNumHeavyAtoms()

    # search in original ligands
    original = 0
    original_sub = 0

    pdb = list(set(pdbs))

    # exact match in original ligand
    if not original_ligands[original_ligands.smiles == smiles].empty:
        original = 1

    # recombined original ligand
    # elif len(pdb) == 1:
    #     print('Original ligand not found:', fragpdbs, pdbs, inchi, smiles)

    # true substructure of original ligands?
    elif not original_ligands[original_ligands.mol >= ligand].empty:
        original_sub = 1

    # exact chembl match
    chembl_match = 0
    chembl_matches = chembl[chembl.smiles == smiles]
    if not chembl_matches.empty:
        chembl_match = 1
    # else:
    #     if original == 1:
    #         print('Original but not ChEMBL:', fragpdbs, pdbs, inchi)

    scaffold = 0
    # Does ligand contain a Hu and Bajorath scaffold as substructure?
    #if not scaffolds[scaffolds.mol <= ligand].empty:
    #    scaffold = 1

    # construct Result object
    result = Result(meta, lipinski, wt, logp, hbd, hba, pains_found, n, original, original_sub, chembl_match, scaffold)

    return result
