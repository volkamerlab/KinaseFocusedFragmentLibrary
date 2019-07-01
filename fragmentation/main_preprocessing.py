from pathlib import Path
from rdkit import Chem

from preprocessing import preprocess_klifs_data, get_folder_name
from discard import is_covalent, contains_phosphate, contains_ribose


path_to_data = Path('../../data/KLIFS_download')
path_to_KLIFS_download = path_to_data / 'overview.csv'
path_to_KLIFS_export = path_to_data / 'KLIFS_export.csv'

KLIFSData = preprocess_klifs_data(path_to_KLIFS_download, path_to_KLIFS_export)
KLIFSData = KLIFSData[KLIFSData.species == 'Human']
before = len(KLIFSData)
KLIFSData = KLIFSData[KLIFSData.dfg == 'in']
after_dfg = len(KLIFSData)

# We are not interested in substrates
KLIFSData = KLIFSData[~KLIFSData.pdb_id.isin(['AMP', 'ADP', 'ATP', 'ACP', 'ANP', 'ADN', 'ADE', 'AGS', 'AN2', 'ANK'])]
after_phosphates = len(KLIFSData)

# count discarded structures
count_ligand_errors = 0
count_pocket_errors = 0
count_multi_ligands = 0
count_riboses = 0
count_covalent = 0
count_dfg_out = before - after_dfg
count_substrates = after_dfg - after_phosphates

filtered_data = KLIFSData.copy()

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = get_folder_name(entry)

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(str(path_to_data / folder / 'ligand.mol2'), removeHs=False)
    pocket = Chem.MolFromMol2File(str(path_to_data / folder / 'pocket.mol2'), removeHs=False)

    try:
        ligandConf = ligand.GetConformer()
    except AttributeError:  # empty molecule
        print('ERROR in ' + folder + ':')
        print('Ligand '+entry.pdb_id+' ('+folder+') could not be loaded. \n')
        count_ligand_errors += 1
        filtered_data = filtered_data.drop(index)
        continue
    try:
        pocketConf = pocket.GetConformer()
    except AttributeError:
        print('ERROR in ' + folder + ':')
        print('Pocket '+folder+' could not be loaded. \n')
        count_pocket_errors += 1
        filtered_data = filtered_data.drop(index)
        continue

    # multiple ligands in one structure
    if '.' in Chem.MolToSmiles(ligand):

        multi_ligands = Chem.GetMolFrags(ligand, asMols=True)
        # do not use phosphate containing ligands
        multi_ligands = [l for l in multi_ligands if not contains_phosphate(l) and not contains_ribose(l)]
        # choose largest one
        sizes = [l.GetNumHeavyAtoms() for l in multi_ligands]
        max_size = max(sizes)

        # if multiple ligands have the same largest size, skip this molecule
        if sizes.count(max_size) > 1:
            print('ERROR in ' + folder + ':')
            print('Ligand consists of multiple molecules. Structure is skipped. \n')
            count_multi_ligands += 1
            filtered_data = filtered_data.drop(index)
            continue

        # if there is a unique largest ligand
        else:
            ligand = multi_ligands[0]
            for l in multi_ligands:
                if l.GetNumHeavyAtoms() > ligand.GetNumHeavyAtoms():
                    ligand = l

    # discard ligands containing phosphates
    if contains_phosphate(ligand):
        print('Phosphate in', entry.pdb, entry.pdb_id, '\n')
        count_substrates += 1
        filtered_data = filtered_data.drop(index)
        continue

    # discard ligands containing riboses
    if contains_ribose(ligand):
        print('Ribose in', entry.pdb, entry.pdb_id)
        count_riboses += 1
        filtered_data = filtered_data.drop(index)
        continue

    # discard covalent ligands
    if is_covalent(entry.pdb, entry.pdb_id, entry.chain):
        print('Covalent inhibitor', entry.pdb, entry.pdb_id, '\n')
        count_covalent += 1
        filtered_data = filtered_data.drop(index)
        continue


filtered_data.to_csv(path_to_data / 'filtered_ligands.csv')


# output statistics
print('\nNumber of discarded structures: ')
print('DFG-out/out-like conformations: ', count_dfg_out)
print('ATP analogs: ', count_substrates)
print('Ribose derivatives: ', count_riboses)
print('Covalent ligands: ', count_covalent)
print('Ligand could not be loaded: ', count_ligand_errors)
print('Pocket could not be loaded: ', count_pocket_errors)
print('Multiple ligands in structure: ', count_multi_ligands)
