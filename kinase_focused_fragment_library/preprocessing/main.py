from pathlib import Path
from rdkit import Chem
import argparse

from preprocessing import preprocess_klifs_data, get_folder_name
from discard import is_covalent, contains_phosphate, contains_ribose, get_ligand_from_multi_ligands

# ============================= INPUT ===============================================

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--klifsdata', type=str, help='path to KLIFS_download folder', required=True)
parser.add_argument('-o', '--fragmentlibrary', type=str, help='output path to fragment library (to write discarded ligands)', required=True)
args = parser.parse_args()

path_to_data = Path(args.klifsdata) / 'KLIFS_download'
path_to_KLIFS_download = path_to_data / 'overview.csv'
path_to_KLIFS_export = path_to_data / 'KLIFS_export.csv'

# ============================= PREPROCESSING ===============================================

# select one structure per PDB
KLIFSData = preprocess_klifs_data(path_to_KLIFS_download, path_to_KLIFS_export)
count_structures = len(KLIFSData)
# select only human kinases
KLIFSData = KLIFSData[KLIFSData.species == 'Human']
before = len(KLIFSData)
# select only DFG-in conformations
KLIFSData = KLIFSData[KLIFSData.dfg == 'in']
after_dfg = len(KLIFSData)
count_dfg_out = before - after_dfg
# We are not interested in Atypical kinases
KLIFSData = KLIFSData[KLIFSData.group != 'Atypical']
count_atypical = after_dfg - len(KLIFSData)

# ============================= INITIALIZATIONS ===============================================

# count discarded structures
count_ligand_errors = 0
count_pocket_errors = 0

errors = []
multi_ligands = []
substrates = []
covalent = []

filtered_data = KLIFSData.copy()

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = get_folder_name(entry)

    # discard substrates
    if entry.pdb_id in ['AMP', 'ADP', 'ATP', 'ACP', 'ANP', 'ADN', 'ADE']:
        filtered_data = filtered_data.drop(index)
        substrates.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
        continue

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
        errors.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
        continue
    try:
        pocketConf = pocket.GetConformer()
    except AttributeError:
        print('ERROR in ' + folder + ':')
        print('Pocket '+folder+' could not be loaded. \n')
        count_pocket_errors += 1
        filtered_data = filtered_data.drop(index)
        errors.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
        continue

    # discard ligands containing phosphates
    if contains_phosphate(ligand):
        print('Phosphate in', entry.pdb, entry.pdb_id, '\n')
        filtered_data = filtered_data.drop(index)
        substrates.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
        continue

    # discard ligands containing riboses
    if contains_ribose(ligand):
        print('Ribose in', entry.pdb, entry.pdb_id)
        filtered_data = filtered_data.drop(index)
        substrates.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
        continue

    # multiple ligands in one structure
    if '.' in Chem.MolToSmiles(ligand):

        ligand = get_ligand_from_multi_ligands(ligand)

        if not ligand:
            print('ERROR in ' + folder + ':')
            print('Ligand consists of multiple molecules. Structure is skipped. \n')
            multi_ligands.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
            continue

    # discard covalent ligands
    if is_covalent(entry.pdb, entry.pdb_id, entry.chain):
        print('Covalent inhibitor', entry.pdb, entry.pdb_id, entry.chain, '\n')
        filtered_data = filtered_data.drop(index)
        covalent.append(entry.pdb+' '+entry.chain+' '+entry.pdb_id)
        continue


filtered_data.to_csv(path_to_data / 'filtered_ligands.csv')


# output statistics
print('Number of starting structures: ', count_structures)
print('\nNumber of discarded structures: ')
print('DFG-out/out-like conformations: ', count_dfg_out)
print('Atypical kinases: ', count_atypical)
print('Ligand could not be loaded: ', count_ligand_errors)
print('Pocket could not be loaded: ', count_pocket_errors)
print('Substrates/substrate analogs: ', len(substrates))
print('Multiple ligands in structure: ', len(multi_ligands))
print('Covalent ligands: ', len(covalent))

folderName = Path(args.fragmentlibrary) / 'discarded_ligands'
if not folderName.exists():
    Path.mkdir(folderName)

with open(folderName / 'errors.txt', 'w') as o:
    for struct in errors:
        o.write(struct+'\n')
with open(folderName / 'multi_ligands.txt', 'w') as o:
    for struct in multi_ligands:
        o.write(struct+'\n')
with open(folderName / 'substrates.txt', 'w') as o:
    for struct in substrates:
        o.write(struct+'\n')
with open(folderName / 'covalent.txt', 'w') as o:
    for struct in covalent:
        o.write(struct+'\n')
