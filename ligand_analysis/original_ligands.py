from rdkit import Chem


def read_original_ligands(frag_dict):

    print('Read original ligands.')

    kinases_pdbs = set()

    for subpocket in frag_dict:

        for frag in frag_dict[subpocket]:
            kinases_pdbs.add((frag.GetProp('kinase'), frag.GetProp('_Name')))

    ligands = []
    for kinase, pdb in kinases_pdbs:
        f = '/home/paula/Masterarbeit/data/KLIFS_download/HUMAN/' + kinase + '/' + pdb + '/ligand.mol2'
        ligand = Chem.MolFromMol2File(f)
        ligand.SetProp('complex_pdb', pdb)
        ligands.append(ligand)

    print('Number of original ligands :', len(ligands))

    return ligands
