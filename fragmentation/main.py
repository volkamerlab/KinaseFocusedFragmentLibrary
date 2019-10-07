from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2
import pandas as pd
import sys
from pathlib import Path
import json

from pocketIdentification import get_subpocket_from_pos, calc_geo_center, fix_small_fragments, calc_subpocket_center
from fragmentation import find_brics_fragments, fragmentation
from classes import Subpocket
from preprocessing import get_folder_name, get_file_name, fix_residue_numbers
from discard import contains_ribose, contains_phosphate
from visualization import visual_subpockets


# ============================= INITIALIZATIONS ===============================================

# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[51], color='0.0, 1.0, 1.0'),  # cyan  # leave out 2? (1m17 will improve)
              Subpocket('AP', residues=[46, 51, 75, 15], color='0.6, 0.1, 0.6'),  # deep purple
              Subpocket('FP', residues=[72, 51, 4, 81], color='0.2, 0.6, 0.2'),  # forest # substituted 4 for 7 and 72 for 74
              Subpocket('GA', residues=[45, 17, 80], color='1.0, 0.5, 0.0'),  # orange
              # Subpocket('BP', residues=[82, 24, 43], color='0.5, 0.0, 1.0')  # purple blue
              Subpocket('B1', residues=[81, 28, 43, 38], color='0.0, 0.5, 1.0'),  # marine
              Subpocket('B2', residues=[18, 24, 70, 83], color='0.5, 0.0, 1.0')  # purple blue # [24, 83, 8, 42] removed 8 because often missing
              ]

# count discarded structures
count_missing_res = 0

count_structures = 0

path_to_library = Path('../FragmentLibrary')

path_to_data = Path('../../data/KLIFS_download')

KLIFSData = pd.read_csv(path_to_data / 'filtered_ligands.csv')

KLIFSData = KLIFSData[KLIFSData.family.isin(['EGFR'])]

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

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = get_folder_name(entry)

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

    # Deal with small fragments
    fix_small_fragments(BRICSFragments, [bond[0] for bond in BRICSBonds])

    # ================================== FRAGMENTATION ==========================================

    skipStructure = False

    # list to store the bonds where we will cleave (as atom tuples)
    bonds = []

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
            bonds.append( ((beginAtom, endAtom), (env_1, env_2)) )

    # print(BRICSFragments[0].environment, Chem.MolToSmiles(BRICSFragments[0].mol))
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

    #     # discard large fragments
    #     if fragment.mol.GetNumHeavyAtoms() > 29:
    #         discardedFragments.append(fragment)
    #     else:
    #         w = Chem.SDWriter(output_files[fragment.subpocket])
    #         w.write(fragment.mol)
    #
    # # ================================ DRAW FRAGMENTS ==========================================
    #
    # # remove Hs and convert to 2D molecules for drawing
    # for fragment in fragments:
    #     fragment.mol = Chem.RemoveHs(fragment.mol)
    #     tmp = AllChem.Compute2DCoords(fragment.mol)
    # img = Draw.MolsToGridImage([fragment.mol for fragment in fragments],
    #                            legends=[fragment.subpocket for fragment in fragments],
    #                            subImgSize=(400, 400))
    # img.save('../output/fragmented_molecules/' + get_file_name(entry) + '.png')

    count_structures += 1


# for subpocket in subpockets:
#     w = Chem.SDWriter(output_files[subpocket.name])
#     output_files[subpocket.name].close()
#     w.close()
#
# # draw discarded fragments
# if discardedFragments:
#     for fragment in discardedFragments:
#         fragment.mol = Chem.RemoveHs(fragment.mol)
#         tmp = AllChem.Compute2DCoords(fragment.mol)
#     img = Draw.MolsToGridImage([fragment.mol for fragment in discardedFragments],
#                                legends=[fragment.structure+' '+fragment.subpocket+' '+str(fragment.mol.GetNumHeavyAtoms())
#                                         for fragment in discardedFragments],
#                                subImgSize=(400, 400), molsPerRow=6)
#     img.save('../output/discarded_fragments.png')
#
# # draw discarded ligands
# if discardedLigands:
#     discardedLigands = [(Chem.RemoveHs(ligand), legend) for ligand, legend in discardedLigands]
#     for ligand in discardedLigands:
#         tmp = AllChem.Compute2DCoords(ligand[0])
#     img = Draw.MolsToGridImage([ligand for ligand, legend in discardedLigands],
#                                legends=[legend for ligand, legend in discardedLigands],
#                                subImgSize=(400, 400), molsPerRow=6)
#     img.save('../output/discarded_ligands.png')


# output statistics
print('Number of fragmented structures: ', count_structures)
print('\nNumber of discarded structures: ')
print('Ligands with too large BRICS fragments: ', len(discardedLigands))
print('Missing residue position could not be inferred: ', count_missing_res)
