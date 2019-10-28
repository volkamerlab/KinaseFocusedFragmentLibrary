from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2
import pandas as pd
import sys
from pathlib import Path
import json

from pocketIdentification import get_subpocket_from_pos, calc_geo_center, fix_small_fragments, calc_subpocket_center, \
    is_valid_subpocket_connection, find_neighboring_fragments
from fragmentation import find_brics_fragments, fragmentation
from classes import Subpocket
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
              # Subpocket('BP', residues=[82, 24, 43], color='0.5, 0.0, 1.0')  # purple blue
              Subpocket('B1', residues=[81, 28, 43, 38], color='0.0, 0.5, 1.0'),  # marine
              Subpocket('B2', residues=[18, 24, 83], color='0.5, 0.0, 1.0')  # purple blue # -70
              ]
se, ap, fp, ga, b1, b2 = subpockets[0], subpockets[1], subpockets[2], subpockets[3], subpockets[4], subpockets[5]

# count discarded structures
count_missing_res = 0

count_structures = 0

path_to_library = Path('../FragmentLibrary')

path_to_data = Path('../../data/KLIFS_download')

KLIFSData = pd.read_csv(path_to_data / 'filtered_ligands.csv')

# KLIFSData = KLIFSData[KLIFSData.kinase.isin(['EGFR', 'AKT1'])]

# clear output files and create output folders
output_files = {}
for subpocket in subpockets:
    folderName = path_to_library / subpocket.name
    if not folderName.exists():
        Path.mkdir(folderName)
    fileName = folderName / (subpocket.name+'.sdf')
    if fileName.exists():
        Path.unlink(fileName)
    output_files[subpocket.name] = fileName.open('a')

discardedFragments = []
discardedLigands = []

invalid_subpocket_connections = {}

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = get_folder_name(entry)

    # special cases where GA can not be disconnected by BRICS which leads to unreasonable connections
    # (but BRICS fragment not large enough to get discarded automatically by the chosen threshold)
    # '3ovv', '3oxt', '3p0m', '3poo': large FP fragment should be in AP, rest is only in FP-II
    if entry.pdb in ['5x5o', '4umt', '4umu', '5mai', '4uyn', '4uzd', '4o0y', '5w5q', '2ycq', '3ovv', '3oxt', '3p0m', '3poo']:
        continue

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(str(path_to_data / folder / 'ligand.mol2'), removeHs=False)
    pocket = Chem.MolFromMol2File(str(path_to_data / folder / 'pocket.mol2'), removeHs=False)

    # multiple ligands in one structure
    if '.' in Chem.MolToSmiles(ligand):

        multi_ligands = Chem.GetMolFrags(ligand, asMols=True)
        # do not use phosphate containing ligands
        multi_ligands = [l for l in multi_ligands if not contains_phosphate(l) and not contains_ribose(l)]
        # choose largest one
        sizes = [l.GetNumHeavyAtoms() for l in multi_ligands]
        max_size = max(sizes)

        # if multiple ligands have the same largest size - should not happen if preprocessing was done correctly
        if sizes.count(max_size) > 1:
            print('ERROR in ' + folder + ':')
            print('Ligand consists of multiple molecules of the same size. \n')
            sys.exit()

        # if there is a unique largest ligand
        else:
            ligand = multi_ligands[0]
            for l in multi_ligands:
                if l.GetNumHeavyAtoms() > ligand.GetNumHeavyAtoms():
                    ligand = l

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

    skipStructure = False

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

    skipStructure = False

    # find BRICS fragments and bonds (as atom numbers)
    BRICSFragments, BRICSBonds = find_brics_fragments(ligand)

    # calculate fragment centers and get nearest subpockets
    for BRICSFragment in BRICSFragments:

        center = calc_geo_center(BRICSFragment.mol.GetAtoms(), BRICSFragment.mol.GetConformer())
        BRICSFragment.center = center

        subpocket = get_subpocket_from_pos(center, subpockets)
        BRICSFragment.subpocket = subpocket

    # discard any ligands where a BRICS fragment is larger than 22 heavy atoms (e.g. staurosporine)
        if BRICSFragment.mol.GetNumHeavyAtoms() > 22:
            discardedLigands.append((ligand, entry.pdb))
            skipStructure = True
            break

    if skipStructure:
        continue

    # ============================= FIX SUBPOCKET ASSIGNMENTS ====================================

    # Adjust subpocket assignments in order to keep small fragments uncleaved
    fix_small_fragments(BRICSFragments, [bond[0] for bond in BRICSBonds])

    # check validity of subpocket connections
    for (beginAtom, endAtom), _ in BRICSBonds:

        firstFragment = next(fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers)
        secondFragment = next(fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers)

        sp_1 = firstFragment.subpocket
        sp_2 = secondFragment.subpocket

        if not is_valid_subpocket_connection(sp_1, sp_2) and sp_1 != sp_2:

            # store invalid subpocket connections
            conn = frozenset((sp_1.name, sp_2.name))
            if conn in invalid_subpocket_connections:
                invalid_subpocket_connections[conn].append(folder)
            else:
                invalid_subpocket_connections[conn] = [folder]

            # # fix invalid subpocket connection
            # # fp-b2 # fp-b1 # ap-b1 # ap-b2 # #se-ga # se-b2
            # # if {sp_1.name, sp_2.name} == {'FP', 'B2'}:
            # #     print(sp_1.name, calc_3d_dist(firstFragment.center, sp_1.center), calc_3d_dist(firstFragment.center, sp_2.center))
            # #     print(sp_2.name, calc_3d_dist(secondFragment.center, sp_2.center), calc_3d_dist(secondFragment.center, sp_1.center))
            # if {sp_1.name, sp_2.name} == {'FP', 'B2'}:
            #
            #     fp_frag = firstFragment if sp_1.name == 'FP' else secondFragment
            #     b2_frag = firstFragment if sp_1.name == 'B2' else secondFragment
            #
            #     # check distance to B2
            #     if calc_3d_dist(b2_frag.center, b2.center) < 5:
            #         print(folder, 'FP-B2 (back pocket)')
            #
            #         # check distance to GA
            #         # print(subpockets[3].name, calc_3d_dist(fp_frag.center, subpockets[3].center), calc_3d_dist(b2_frag.center, subpockets[3].center))
            #         if calc_3d_dist(fp_frag.center, ga.center) < 7:
            #             # fp_frag.subpocket = ga
            #             # # print('FP -> GA')
            #             # # check size of new GA fragment
            #             # if fp_frag.mol.GetNumHeavyAtoms() < 3:
            #             #     print(folder, 'fp-b2, ga too small')
            #             print(folder, 'FP-GA-B2: FP -> GA')
            #
            #         elif calc_3d_dist(b2_frag.center, ga.center) < 7:
            #             # b2_frag.subpocket = ga
            #             # # print('B2 -> GA')
            #             # if b2_frag.mol.GetNumHeavyAtoms() < 3:
            #             #     print(folder, 'fp-b2, ga too small')
            #             print(folder, 'FP-GA-B2: B2 -> GA')
            #
            #     # if B2 fragment is not close enough to B2, put all in FP
            #     else:
            #         print(folder, 'FP-B2: B2 -> FP')
            #         # assign B2 fragments to FP
            #         b2_frag.subpocket = fp
            #         num_fixed_fragments = 1
            #         fixed_frag = b2_frag
            #         while num_fixed_fragments > 0:
            #             num_fixed_fragments = 0
            #             for nf in find_neighboring_fragments(fixed_frag, BRICSFragments, [bond[0] for bond in BRICSBonds]):
            #                 if nf.subpocket == b2:
            #                     nf.subpocket = fp
            #                     num_fixed_fragments += 1
            #                     fixed_frag = nf

        # # check special cases SE-FP and FP-GA
        # elif {sp_1.name, sp_2.name} == {'SE', 'FP'} or {sp_1.name, sp_2.name} == {'FP', 'GA'}:
        #
        #     conn = frozenset((sp_1.name, sp_2.name))
        #     if conn in invalid_subpocket_connections:
        #         invalid_subpocket_connections[conn].append(folder)
        #     else:
        #         invalid_subpocket_connections[conn] = [folder]

    # ================================== FRAGMENTATION ==========================================

    skipStructure = False

    # list to store the bonds where we will cleave (as atom tuples)
    bonds = []
    count = 0
    # iterate over BRICS bonds
    for (beginAtom, endAtom), (env_1, env_2) in BRICSBonds:

        # find corresponding fragments
        firstFragment = next(fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers)
        secondFragment = next(fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers)

        # get environment types of the brics fragments
        firstFragment.environment = env_1
        secondFragment.environment = env_2

        # check if subpockets differ
        if firstFragment.subpocket != secondFragment.subpocket:
            # store this bond as a bond where we will cleave
            bonds.append((beginAtom, endAtom))

    # check validity of subpocket connections
    for (beginAtom, endAtom) in bonds:

        firstFragment = next(fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers)
        secondFragment = next(fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers)

        sp_1 = firstFragment.subpocket
        sp_2 = secondFragment.subpocket

        if not is_valid_subpocket_connection(sp_1, sp_2):

            print(folder, 'Invalid subpocket connection:', sp_1.name, '-', sp_2.name)

        # discard this ligand if SE and GA subpockets are connected
        if {sp_1.name, sp_2.name} == {'SE', 'GA'}:

            skipStructure = True
            break

    if skipStructure:
        continue

    # actual fragmentation
    fragments = fragmentation(ligand, bonds, BRICSFragments)

    # ================================ FRAGMENT LIBRARY ========================================

    # add fragments to their respective pool

    for fragment in fragments:
        # store PDB where this fragment came from
        fragment.structure = get_file_name(entry)

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

        # discard large fragments
        if fragment.mol.GetNumHeavyAtoms() > 29:
            discardedFragments.append(fragment)
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


for subpocket in subpockets:
    w = Chem.SDWriter(output_files[subpocket.name])
    output_files[subpocket.name].close()
    w.close()

# draw discarded fragments
if discardedFragments:
    for fragment in discardedFragments:
        fragment.mol = Chem.RemoveHs(fragment.mol)
        tmp = AllChem.Compute2DCoords(fragment.mol)
    img = Draw.MolsToGridImage([fragment.mol for fragment in discardedFragments],
                               legends=[fragment.structure+' '+fragment.subpocket.name+' '+str(fragment.mol.GetNumHeavyAtoms())
                                        for fragment in discardedFragments],
                               subImgSize=(400, 400), molsPerRow=6)
    img.save('../output/discarded_fragments.png')

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

# print invalid subpocket connections
for conn in invalid_subpocket_connections:
    print([sp for sp in conn], len(invalid_subpocket_connections[conn]), len(invalid_subpocket_connections[conn])/count_structures*100)
    for struct in invalid_subpocket_connections[conn]:
        print(struct)

# print('\nSE-FP-GA')
# for struct in invalid_subpocket_connections[frozenset(('SE', 'FP'))]:
#     for struct2 in invalid_subpocket_connections[frozenset(('FP', 'GA'))]:
#         if struct == struct2:
#             print(struct)
