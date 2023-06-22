from pathlib import Path
from rdkit import Chem
import argparse
import pandas as pd

from kinase_focused_fragment_library.preprocessing.preprocessing import \
    read_klifs_meta_data, choose_best_klifs_structure, get_folder_name
from kinase_focused_fragment_library.preprocessing.discard import \
    is_covalent, contains_phosphate, contains_ribose, get_ligand_from_multi_ligands


def main():

    # ============================= INPUT ===============================================

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--klifsdata', type=str, help='path to KLIFS_download folder', required=True)
    parser.add_argument('-o', '--fragmentlibrary', type=str, help='output path to fragment library (to write discarded ligands)', required=True)
    args = parser.parse_args()

    path_to_data = Path(args.klifsdata)
    path_to_KLIFS_download = path_to_data / 'data.csv'
    path_to_KLIFS_export = path_to_data / 'data.csv'

    # ============================= READ METADATA ===============================================

    # select one structure per PDB
    KLIFSData = read_klifs_meta_data(path_to_KLIFS_download, path_to_KLIFS_export)
    count_structures = len(KLIFSData)

    # ================== Filter by kinase group, species, DFG-conformation ======================

    # We are not interested in Atypical kinases
    # KLIFSData = KLIFSData[KLIFSData.group != 'Atypical']
    count_atypical = count_structures - len(KLIFSData)
    # select only human kinases
    KLIFSData = KLIFSData[KLIFSData.species == 'Human']
    before = len(KLIFSData)
    # select only DFG-in conformations
    # KLIFSData = KLIFSData[KLIFSData.dfg == 'in']
    after_dfg = len(KLIFSData)
    count_dfg_out = before - after_dfg

    # ======================= SELECT ONE STRUCTURE PER PDB =======================================

    # select best structure per PDB
    # KLIFSData = choose_best_klifs_structure(KLIFSData)
    count_structures = len(KLIFSData)
    count_duplicate_pdbs = after_dfg - count_structures

    # ============================= INITIALIZATIONS ==============================================

    # count discarded structures
    count_ligand_errors = 0
    count_pocket_errors = 0
    count_multi_ligands = 0
    count_substrates = 0
    count_covalent = 0

    # output file with metadata of structures chosen for fragmentation
    filtered_data = KLIFSData.copy()

    # output file with metadata of discarded structures
    discarded_structures = []

    # =========================== ITERATE OVER STRUCTURES =========================================

    # iterate over molecules
    for index, entry in KLIFSData.iterrows():

        folder = get_folder_name(entry)

        # ============================ FILTER SUBSTRATES ==========================================

        # discard substrates
        if entry.pdb in ['AMP', 'ADP', 'ATP', 'ACP', 'ANP', 'ADN', 'ADE']:
            filtered_data = filtered_data.drop(index)
            count_substrates += 1
            entry['violation'] = 'Substrate'
            discarded_structures.append(entry)
            continue

        # ==================== FILTER UNLOADBALE STRUCTURES =======================================

        # load ligand and binding pocket to rdkit molecules
        try:
            ligand_path = str(path_to_data / 'ligand_poses' / folder / 'ligand.mol2')
            ligand = Chem.MolFromMol2File(ligand_path, removeHs=False, sanitize=False)
            pocket = Chem.MolFromMol2File(str(path_to_data / 'ligand_poses' / folder / 'pocket.mol2'), removeHs=False)
        except OSError:
            print('ERROR in ' + folder + ':')
            print('Ligand or pocket file '+str(entry.ident) + ' from ', ligand_path)#+' ('+folder+') could not be loaded. \n')
            count_ligand_errors += 1
            filtered_data = filtered_data.drop(index)
            entry['violation'] = 'Unloadable ligand or pocket file'
            discarded_structures.append(entry)
            continue

        # check if KLIFS ligand and protein are loadable with RDKit
        try:
            ligandConf = ligand.GetConformer()
        except AttributeError:  # empty molecule
            print('ERROR in ' + folder + ':')
            print('Ligand '+entry.pdb+' ('+folder+') could not be loaded. \n')
            count_ligand_errors += 1
            filtered_data = filtered_data.drop(index)
            entry['violation'] = 'Unloadable ligand'
            discarded_structures.append(entry)
            continue
        try:
            pocketConf = pocket.GetConformer()
        except AttributeError:
            print('ERROR in ' + folder + ':')
            print('Pocket '+folder+' could not be loaded. \n')
            count_pocket_errors += 1
            filtered_data = filtered_data.drop(index)
            entry['violation'] = 'Unloadable pocket'
            discarded_structures.append(entry)
            continue

        # ===================== FILTER PHOSPHATES AND RIBOSES ===================================

        # discard ligands containing phosphates
        if contains_phosphate(ligand):
            print('Phosphate in', entry.pdb, entry.pdb, '\n')
            filtered_data = filtered_data.drop(index)
            count_substrates += 1
            entry['violation'] = 'Contains phosphate'
            discarded_structures.append(entry)
            continue

        # discard ligands containing riboses
        if contains_ribose(ligand):
            print('Ribose in', entry.pdb, entry.pdb)
            filtered_data = filtered_data.drop(index)
            count_substrates += 1
            entry['violation'] = 'Contains ribose'
            discarded_structures.append(entry)
            continue

        # =================== FILTER MULTIPLE LIGANDS PER STRUCTURE ==============================

        # multiple ligands in one structure
        if '.' in Chem.MolToSmiles(ligand):

            ligand = get_ligand_from_multi_ligands(ligand)

            if not ligand:
                print('ERROR in ' + folder + ':')
                print('Ligand consists of multiple molecules. Structure is skipped. \n')
                filtered_data = filtered_data.drop(index)
                count_multi_ligands += 1
                entry['violation'] = 'Multiple ligands present'
                discarded_structures.append(entry)
                continue

        # ============================ FILTER COVALENT LIGANDS ===================================

        # discard covalent ligands
        if is_covalent(entry.ligand_pdb, entry.pdb, entry.chain):
            print('Covalent inhibitor', entry.ligand_pdb, entry.pdb, entry.chain, '\n')
            filtered_data = filtered_data.drop(index)
            count_covalent += 1
            entry['violation'] = 'Covalent ligand'
            discarded_structures.append(entry)
            continue


    # ============================= OUTPUT ===============================================

    # write to output files
    filtered_data.to_csv(path_to_data / 'filtered_ligands.csv')

    folderName = Path(args.fragmentlibrary) / 'discarded_ligands'
    Path.mkdir(folderName, parents=True, exist_ok=True)
    discarded_structures = pd.DataFrame(discarded_structures).reset_index(drop=True)
    discarded_structures.to_csv(folderName / 'preprocessing.csv')

    # output statistics
    print('Number of starting structures: ', count_structures)
    print('\nNumber of discarded structures: ', len(discarded_structures))
    print('DFG-out/out-like conformations: ', count_dfg_out)
    print('Atypical kinases: ', count_atypical)
    print('Ligand could not be loaded: ', count_ligand_errors)
    print('Pocket could not be loaded: ', count_pocket_errors)
    print('Substrates/substrate analogs: ', count_substrates)
    print('Multiple ligands in structure: ', count_multi_ligands)
    print('Covalent ligands: ', count_covalent)


if __name__ == "__main__":
    main()
