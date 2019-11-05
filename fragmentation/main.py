from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2
import pandas as pd
import sys
from pathlib import Path
import json

from pocketIdentification import get_subpocket_from_pos, calc_geo_center, fix_small_fragments, calc_subpocket_center, \
    is_valid_subpocket_connection
from fragmentation import find_brics_fragments, fragment_between_atoms, set_atom_properties
from classes import Subpocket, Fragment
from preprocessing import get_folder_name, get_file_name, fix_residue_numbers
from discard import contains_ribose, contains_phosphate
from visualization import visual_subpockets
from functions import calc_3d_dist


# ============================= INITIALIZATIONS ===============================================

# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[51], color='0.0, 1.0, 1.0'),  # cyan  # leave out 2? (1m17 will improve)
              Subpocket('AP', residues=[46, 51, 75, 15], color='0.6, 0.1, 0.6'),  # deep purple
              Subpocket('FP', residues=[72, 51, 10, 81], color='0.2, 0.6, 0.2'),  # forest # 4 -> 10
              Subpocket('GA', residues=[45, 17, 81], color='1.0, 0.5, 0.0'),  # orange  # 80 -> 82 -> 81
              Subpocket('B1', residues=[81, 28, 43, 38], color='0.0, 0.5, 1.0'),  # marine
              Subpocket('B2', residues=[18, 24, 70, 83], color='0.5, 0.0, 1.0'),  # purple blue # -70
              ]
se, ap, fp, ga, b1, b2 = subpockets[0], subpockets[1], subpockets[2], subpockets[3], subpockets[4], subpockets[5]

# count discarded structures
count_missing_res = 0
count_not_ap = 0
count_structures = 0

path_to_library = Path('../FragmentLibrary')

path_to_data = Path('../../data/KLIFS_download')

KLIFSData = pd.read_csv(path_to_data / 'filtered_ligands.csv')

# KLIFSData = KLIFSData[KLIFSData.kinase.isin(['TTK'])]

# clear output files and create output folders
output_files = {}
for subpocket in subpockets+[Subpocket('X')]:
    folderName = path_to_library / subpocket.name
    if not folderName.exists():
        Path.mkdir(folderName)
    fileName = folderName / (subpocket.name+'.sdf')
    if fileName.exists():
        Path.unlink(fileName)
    output_files[subpocket.name] = fileName.open('a')

discardedLigands = []

invalid_subpocket_connections = {}
distances = {}
outliers = set()

fragment_sizes = dict()
for size in range(1, 40):
    fragment_sizes[size] = 0

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = get_folder_name(entry)

    skipStructure = False

    # special cases where GA or FP can not be disconnected by BRICS which leads to unreasonable subpocket assignments
    # if entry.pdb in ['5mai', '4o0y', '5w5q']:
    #     continue

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(str(path_to_data / folder / 'ligand.mol2'), removeHs=False)
    pocket = Chem.MolFromMol2File(str(path_to_data / folder / 'pocket.mol2'), removeHs=False)

    # multiple ligands in one structure
    if '.' in Chem.MolToSmiles(ligand):

        multi_ligands = Chem.GetMolFrags(ligand, asMols=True)
        # do not use structures including substrates
        phosphate_ligands = [l for l in multi_ligands if contains_phosphate(l) or contains_ribose(l)]
        if phosphate_ligands:
            print('ERROR in ' + folder + ':')
            print('Ligand consists of multiple molecules. Structure is skipped. \n')
            continue
        # get only large ligands
        multi_ligands = [l for l in multi_ligands if l.GetNumHeavyAtoms() > 14]
        # if there is more than one large ligand, discard this structure
        if len(multi_ligands) != 1:
            print('ERROR in ' + folder + ':')
            print('Ligand consists of multiple molecules. Structure is skipped. \n')
            continue
        else:
            ligand = multi_ligands[0]

    lenLigand = ligand.GetNumAtoms()

    # read atom information from binding pocket mol2 file (necessary for residue information)
    pocketMol2 = PandasMol2().read_mol2(str(path_to_data / folder / 'pocket.mol2'),
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary_structure', str)}).df

    # convert string to list
    missing_residues = json.loads(entry.missing_residues)
    # fix residue IDs
    pocketMol2 = fix_residue_numbers(pocketMol2, missing_residues)

    # ============================ SUBPOCKET CENTERS =========================================

    # calculate subpocket centers
    for subpocket in subpockets:

        subpocket.center = calc_subpocket_center(subpocket, pocket, pocketMol2, folder)
        # skip structure if no center could be calculated because of missing residues
        if subpocket.center is None:
            count_missing_res += 1
            skipStructure = True
            break

    # skip this molecule if important residues are missing
    if skipStructure:
        continue

    # visualize subpocket centers using PyMOL
    visual_subpockets(subpockets, folder)

    # ================================ BRICS FRAGMENTS ==========================================

    # find BRICS fragments and bonds (as atom numbers)
    BRICSFragments, BRICSBonds = find_brics_fragments(ligand)

    # calculate fragment centers and get nearest subpockets
    for BRICSFragment in BRICSFragments:

        center = calc_geo_center(BRICSFragment.mol.GetAtoms(), BRICSFragment.mol.GetConformer())
        BRICSFragment.center = center

        subpocket, distance = get_subpocket_from_pos(center, subpockets, distances)

        # check distance to nearest subpocket
        if distance >= 8:
            # if distance to nearest subpocket is too large, assign this fragment to X ("bin" pocket)
            BRICSFragment.subpocket = Subpocket('X-'+subpocket.name)
        else:
            # else assign this fragment to its nearest subpocket
            BRICSFragment.subpocket = subpocket

    # discard any ligands where a BRICS fragment is larger than 22 heavy atoms (e.g. staurosporine)
        if BRICSFragment.mol.GetNumHeavyAtoms() > 22:
            discardedLigands.append((ligand, entry.pdb))
            skipStructure = True
            break

    if skipStructure:
        continue

    # Adjust subpocket assignments in order to keep small fragments uncleaved
    fix_small_fragments(BRICSFragments, [bond[0] for bond in BRICSBonds], 4)

    # ================================== FRAGMENTATION ==========================================

    # list to store the bonds where we will cleave
    atom_tuples = []
    count = 0
    # iterate over BRICS bonds
    for (beginAtom, endAtom), (env_1, env_2) in BRICSBonds:

        # find corresponding fragments
        firstFragment = next(fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers)
        secondFragment = next(fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers)

        # set environment types of the brics fragments
        firstFragment.environment = env_1
        secondFragment.environment = env_2

        # check if subpockets differ
        if firstFragment.subpocket != secondFragment.subpocket:
            # store this bond as a bond where we will cleave
            atom_tuples.append((beginAtom, endAtom))

    # fragmentation of the ligand at the calculated bonds that separate two subpockets
    fragment_mols, fragment_atoms = fragment_between_atoms(ligand, atom_tuples)

    # ============================= SUBPOCKET ASSIGNMENTS ====================================

    # iterate over new fragments and create Fragment objects
    fragments = []
    for (atomNumbers, mol) in zip(fragment_atoms, fragment_mols):

        # get subpocket corresponding to fragment (Is there a better way?)
        subpocket = next(brics_fragment.subpocket for brics_fragment in BRICSFragments if atomNumbers[0] in brics_fragment.atomNumbers)
        # create Fragment object
        fragments.append(Fragment(mol=mol, atomNumbers=atomNumbers, subpocket=subpocket))

    # skip this structure if it does not contain an AP fragment
    if ap not in [fragment.subpocket for fragment in fragments]:
        count_not_ap += 1
        continue

    # check for FP-B2 connections
    for (beginAtom, endAtom) in atom_tuples:

        firstFragment = next(fragment for fragment in fragments if beginAtom in fragment.atomNumbers)
        secondFragment = next(fragment for fragment in fragments if endAtom in fragment.atomNumbers)

        sp_1 = firstFragment.subpocket
        sp_2 = secondFragment.subpocket

        # if FP and B2 are connected, check distance to GA and if close enough put FP fragment to GA
        if {sp_1.name, sp_2.name} == {'FP', 'B2'}:

            fp_frag = firstFragment if sp_1.name == 'FP' else secondFragment
            b2_frag = firstFragment if sp_1.name == 'B2' else secondFragment
            fp_frag.center = calc_geo_center(fp_frag.mol.GetAtoms(), fp_frag.mol.GetConformer())

            # if FP fragment is close to GA
            if calc_3d_dist(fp_frag.center, ga.center) < 5:
                # assign FP fragment to GA pocket
                fp_frag.subpocket = ga

            # if FP is not close to GA, there is an actual FP-B2 connection which we do not want to have
            # hence we assign the B2 fragment to X
            else:
                b2_frag.subpocket = Subpocket('X-B2')

    # check subpocket connections again after adaptation
    for (beginAtom, endAtom) in atom_tuples:

        firstFragment = next(fragment for fragment in fragments if beginAtom in fragment.atomNumbers)
        secondFragment = next(fragment for fragment in fragments if endAtom in fragment.atomNumbers)

        sp_1 = firstFragment.subpocket
        sp_2 = secondFragment.subpocket

        if not is_valid_subpocket_connection(sp_1, sp_2):

            print(folder+':')
            print('Invalid subpocket connection:', sp_1.name, '-', sp_2.name, '. Structure is skipped.\n')

            # store invalid subpocket connections
            conn = frozenset((sp_1.name, sp_2.name))
            if conn in invalid_subpocket_connections:
                invalid_subpocket_connections[conn].append(folder)
            else:
                invalid_subpocket_connections[conn] = [folder]

            skipStructure = True
            break

    if skipStructure:
        continue

    # set atom properties: atom ids, subpockets, and BRICS environments
    set_atom_properties(fragments, atom_tuples, BRICSFragments)

    # ================================ FRAGMENT LIBRARY ========================================

    # add fragments to their respective pool

    for fragment in fragments:

        # store PDB where this fragment came from
        fragment.structure = get_file_name(entry)

        fragment_sizes[fragment.mol.GetNumHeavyAtoms()] += 1

        # fragment properties
        # this sets the PDB code as 'name' of the fragment at the top of the SD file entry
        fragment.mol.SetProp('_Name', fragment.structure)
        # set other properties to be stored in the SD file
        fragment.mol.SetProp('kinase', entry.kinase)
        fragment.mol.SetProp('family', entry.family)
        fragment.mol.SetProp('group', entry.group)
        fragment.mol.SetProp('complex_pdb', entry.pdb)
        fragment.mol.SetProp('ligand_pdb', entry.pdb_id)
        fragment.mol.SetProp('alt', entry.alt)
        fragment.mol.SetProp('chain', entry.chain)

        # atom properties as fragment property
        Chem.CreateAtomStringPropertyList(fragment.mol, 'subpocket')
        Chem.CreateAtomStringPropertyList(fragment.mol, 'environment')

        if fragment.subpocket.name.startswith('X'):
            w = Chem.SDWriter(output_files['X'])
            outliers.add(folder)
        else:
            w = Chem.SDWriter(output_files[fragment.subpocket.name])
        w.write(fragment.mol)

    # ================================ DRAW FRAGMENTS ==========================================

    # remove Hs and convert to 2D molecules for drawing
    for fragment in fragments:
        fragment.mol = Chem.RemoveHs(fragment.mol)
        tmp = AllChem.Compute2DCoords(fragment.mol)

    img = Draw.MolsToGridImage([fragment.mol for fragment in fragments],
                               legends=[fragment.subpocket.name for fragment in fragments],
                               subImgSize=(400, 400))
    img.save('../output/fragmented_molecules/' + get_file_name(entry) + '.png')

    count_structures += 1


for subpocket in subpockets+[Subpocket('X')]:
    w = Chem.SDWriter(output_files[subpocket.name])
    output_files[subpocket.name].close()
    w.close()

# draw discarded ligands
if discardedLigands:
    discardedLigands = [(Chem.RemoveHs(ligand), legend) for ligand, legend in discardedLigands]
    for ligand in discardedLigands:
        tmp = AllChem.Compute2DCoords(ligand[0])
    img = Draw.MolsToGridImage([ligand for ligand, legend in discardedLigands],
                               legends=[legend for ligand, legend in discardedLigands],
                               subImgSize=(400, 400), molsPerRow=6)
    img.save('../output/discarded_ligands.png')


# output statistics
print('Number of fragmented structures: ', count_structures)
print('\nNumber of discarded structures: ')
print('Ligands with too large BRICS fragments: ', len(discardedLigands))
print('Missing residue position could not be inferred: ', count_missing_res)
print('Ligands not occupying AP:', count_not_ap)

# print invalid subpocket connections
for conn in invalid_subpocket_connections:
    print([sp for sp in conn], len(invalid_subpocket_connections[conn]), len(invalid_subpocket_connections[conn])/count_structures*100)
    for struct in invalid_subpocket_connections[conn]:
        print(struct)

print('Fragment sizes:')
print(fragment_sizes)
