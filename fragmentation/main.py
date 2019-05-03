from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2

from pocketIdentification import get_subpocket_from_pos, calc_geo_center, fix_small_fragments, calc_subpocket_center
from fragmentation import find_brics_fragments, fragmentation
from classes import Subpocket
from preprocessing import preprocess_klifs_data, get_folder_name, get_file_name, fix_residue_numbers
from visualization import visual_subpockets

from pathlib import Path

# ============================= INITIALIZATIONS ===============================================

# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[51], color='0.0, 1.0, 1.0'),  # cyan  # leave out 2? (1m17 will improve)
              Subpocket('AP', residues=[46, 51, 75, 15], color='0.6, 0.1, 0.6'),  # deeppurple
              Subpocket('FP', residues=[72, 51, 4, 81], color='0.2, 0.6, 0.2'),  # forest # substituted 4 for 7 and 72 for 74
              Subpocket('GA', residues=[45, 17, 80], color='1.0, 0.5, 0.0'),  # orange
              # Subpocket('BP', residues=[82, 24, 43], color='0.5, 0.0, 1.0')  # purple blue
              Subpocket('B1', residues=[81, 28, 43, 38], color='0.0, 0.5, 1.0'),  # marine
              Subpocket('B2', residues=[18, 24, 70, 83], color='0.5, 0.0, 1.0')  # purple blue # [24, 83, 8, 42] removed 8 because often missing
              ]

# count discarded structures
count_ligand_errors = 0
count_pocket_errors = 0
count_multi_ligands = 0
count_missing_res = 0
# count_phosphates = 0
# count_dfg_out = 0
count_structures = 0

# ============================= DATA PREPARATION ============================================

path_to_library = Path('../FragmentLibrary')

path_to_data = Path('../../data/KLIFS_download')
path_to_KLIFS_download = path_to_data / 'overview.csv'
path_to_KLIFS_export = path_to_data / 'KLIFS_export.csv'

KLIFSData = preprocess_klifs_data(path_to_KLIFS_download, path_to_KLIFS_export)
KLIFSData = KLIFSData[KLIFSData.species == 'Human']
before = len(KLIFSData)
KLIFSData = KLIFSData[KLIFSData.dfg == 'in']
after_dfg = len(KLIFSData)
count_dfg_out = before - after_dfg
# We are not interested in adenosine phosphates
KLIFSData = KLIFSData[~KLIFSData.pdb_id.isin(['AMP', 'ADP', 'ATP', 'ACP', 'ANP'])]
after_phosphates = len(KLIFSData)
count_phosphates = after_dfg - after_phosphates

# KLIFSData = KLIFSData[KLIFSData.family.isin(['RAF', 'EGFR', 'CDK'])]


# clear output files and create output folders
for subpocket in subpockets:
    folderName = path_to_library / subpocket.name
    if not folderName.exists():
        Path.mkdir(folderName)
    kinases = set(KLIFSData.kinase)
    for kinase in kinases:
        fileName = folderName / (subpocket.name+'.sdf')
        if fileName.exists():
            Path.unlink(fileName)

discardedFragments = []
discardedLigands = []

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = get_folder_name(entry)
    # print(folder, entry.dfg, entry.ac_helix)

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(str(path_to_data / folder / 'ligand.mol2'), removeHs=False)
    pocket = Chem.MolFromMol2File(str(path_to_data / folder / 'pocket.mol2'), removeHs=False)

    try:
        ligandConf = ligand.GetConformer()
    except AttributeError:  # empty molecule
        print('ERROR in ' + folder + ':')
        print('Ligand '+entry.pdb_id+' ('+folder+') could not be loaded. \n')
        count_ligand_errors += 1
        continue
    try:
        pocketConf = pocket.GetConformer()
    except AttributeError:
        print('ERROR in ' + folder + ':')
        print('Pocket '+folder+' could not be loaded. \n')
        count_pocket_errors += 1
        continue

    # multiple ligands in one structure
    if '.' in Chem.MolToSmiles(ligand):
        # choose largest one
        multi_ligands = Chem.GetMolFrags(ligand, asMols=True)
        sizes = [l.GetNumHeavyAtoms() for l in multi_ligands]
        max_size = max(sizes)

        # if multiple ligands have the same largest size, skip this molecule
        if sizes.count(max_size) > 1:
            print('ERROR in ' + folder + ':')
            print('Ligand consists of multiple molecules. Structure is skipped. \n')
            count_multi_ligands += 1
            continue

        # if there is a unique largest ligand
        else:
            ligand = multi_ligands[0]
            for l in multi_ligands:
                if l.GetNumHeavyAtoms() > ligand.GetNumHeavyAtoms():
                    ligand = l

    lenLigand = ligand.GetNumAtoms()

    # read atom information from binding pocket mol2 file (necessary for residue information)
    # pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')
    pocketMol2 = PandasMol2().read_mol2(str(path_to_data / folder / 'pocket.mol2'),
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary_structure', str)}).df
    # fix residue IDs
    pocketMol2 = fix_residue_numbers(pocketMol2, entry.missing_residues)

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

    # # get subpocket for each ligand atom
    # for a, atom in enumerate(ligand.GetAtoms()):
    #     # subpocket = getSubpocketFromAtomDistances(a, ligandConf, pocketConf, residues)
    #     subpocket = getSubpocketFromAtom(a, ligandConf, subpockets)
    #     atom.SetProp('subpocket', subpocket)

    # ================================ BRICS FRAGMENTS ==========================================

    skipStructure = False

    # find BRICS fragments and bonds (as atom numbers)
    BRICSFragments, BRICSBonds = find_brics_fragments(ligand)

    # calculate fragment centers and get nearest subpockets
    for BRICSFragment in BRICSFragments:

        center = calc_geo_center(BRICSFragment.mol.GetAtoms(), BRICSFragment.mol.GetConformer())
        BRICSFragment.center = center

    # ========================== SUBPOCKET IDENTIFICATION ========================================
    #
    #     # ---------------------------------------
    #     # instead of getSubpocketFromPos() function in order to find ambiguous fragments
    #
    #     # calculate distances from subpockets to fragments
    #     smallestDistance = sys.maxsize  # set smallest distance as max integer value
    #     nearestSubpocket = Subpocket('noSubpocket')
    #     distances = []
    #     for subpocket in subpockets:
    #         distance = calculate3DDistance(center, subpocket.center)
    #         distances.append(distance)
    #         if distance < smallestDistance:
    #             nearestSubpocket = subpocket
    #             smallestDistance = distance
    #
    #     # draw ambiguous fragments
    #     minDist = min(distances)
    #     for d, dist in enumerate(distances):
    #         # if there is a value near the minimum distance
    #         if minDist - 0.7 < dist < minDist + 0.7 and subpockets[d] != nearestSubpocket:
    #             label = subpockets[d].name+'+'+nearestSubpocket.name
    #             Draw.MolToFile(BRICSFragment.mol, '../ambiguous_fragments/'+getFileName(entry)+'_'+label+'.png')
    #
    #     BRICSFragment.subpocket = nearestSubpocket.name
    #
    #     # -------------------------------------------

        subpocket = get_subpocket_from_pos(center, subpockets)
        BRICSFragment.subpocket = subpocket

    # discard any ligands where a BRICS fragment is larger than 22 heavy atoms (e.g. staurosporine)
        if BRICSFragment.mol.GetNumHeavyAtoms() > 22:
            discardedLigands.append([ligand, entry.pdb])
            skipStructure = True
            break

    if skipStructure:
        continue

    # Deal with small fragments
    fix_small_fragments(BRICSFragments, BRICSBonds)

    # ================================== FRAGMENTATION ==========================================

    skipStructure = False

    # list to store the bonds where we will cleave (as atom tuples)
    bonds = []

    # iterate over BRICS bonds
    for beginAtom, endAtom in BRICSBonds:

        # find corresponding fragments
        firstFragment = [fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers][0]
        secondFragment = [fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers][0]

        # # check validity of subpockets
        # if not checkSubpockets(firstFragment.subpocket, secondFragment.subpocket):
        #
        #     print('ERROR in '+folder+':')
        #     print("Subpockets "+firstFragment.subpocket+" and "+secondFragment.subpocket+" can not be connected."
        #                                                                                  "Structure is skipped. \n")
        #     # skip this molecule if subpocket definition is not valid
        #     skipStructure = True
        #     break

        # check if subpockets differ
        if firstFragment.subpocket != secondFragment.subpocket:
            # store this bond as a bond where we will cleave
            bonds.append((beginAtom, endAtom))

    # # skip this molecule if subpocket definition is not valid
    # if skipStructure:
    #     continue

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
        fragment.mol.SetProp('group', entry.groups)
        fragment.mol.SetProp('complex_pdb', entry.pdb)
        fragment.mol.SetProp('ligand_pdb', entry.pdb_id)
        fragment.mol.SetProp('alt', entry.alt)
        fragment.mol.SetProp('chain', entry.chain)

        # atom properties as fragment property
        # RemoveHs() does not remove all hydrogens, inconsistency with removeHs=True flag
        # fragment.mol = Chem.RemoveHs(fragment.mol)  # remove hydrogens for consistency reasons (when using PandasTools.LoadSDF)
        # create list of atom properties
        # Chem.CreateAtomStringPropertyList(fragment.mol, 'neighboringSubpocket')
        Chem.CreateAtomStringPropertyList(fragment.mol, 'subpocket')

        # discard large fragments
        if fragment.mol.GetNumHeavyAtoms() > 29:
            discardedFragments.append(fragment)
        else:
            # output_file = path_to_library+fragment.subpocket+'/'+getFileName(entry)+'.sdf'
            output_file = (path_to_library / fragment.subpocket / (fragment.subpocket+'.sdf')).open('a')
            # print(Chem.MolToMolBlock(fragment.mol), file=open(output_file, 'a'))
            w = Chem.SDWriter(output_file)
            w.write(fragment.mol)

    # ================================ DRAW FRAGMENTS ==========================================

    # convert to 2D molecules for drawing
    for fragment in fragments:
        tmp = AllChem.Compute2DCoords(fragment.mol)

    img = Draw.MolsToGridImage([Chem.RemoveHs(fragment.mol) for fragment in fragments],
                               legends=[fragment.subpocket for fragment in fragments],
                               subImgSize=(400, 400))
    img.save('../output/fragmented_molecules/' + get_file_name(entry) + '.png')

    count_structures += 1

# draw discarded fragments
if discardedFragments:
    img = Draw.MolsToGridImage([Chem.RemoveHs(fragment.mol) for fragment in discardedFragments],
                               legends=[fragment.structure+' '+fragment.subpocket+' '+str(fragment.mol.GetNumHeavyAtoms())
                                        for fragment in discardedFragments],
                               subImgSize=(400, 400), molsPerRow=6)
    img.save('../output/discarded_fragments.png')

# draw discarded ligands
if discardedLigands:
    for ligand in discardedLigands:
        tmp = AllChem.Compute2DCoords(ligand[0])
    img = Draw.MolsToGridImage([Chem.RemoveHs(ligand[0]) for ligand in discardedLigands],
                               legends=[ligand[1] for ligand in discardedLigands],
                               subImgSize=(400, 400), molsPerRow=6)
    img.save('../output/discarded_ligands.png')


# output statistics
print('Number of fragmented structures: ', count_structures)
print('\nNumber of discarded structures: ')
print('DFG-out/out-like conformations: ', count_dfg_out)
print('A*P ligands: ', count_phosphates)
print('Ligand could not be loaded: ', count_ligand_errors)
print('Pocket could not be loaded: ', count_pocket_errors)
print('Multiple ligands in structure: ', count_multi_ligands)
print('Ligands with too large BRICS fragments: ', len(discardedLigands))
print('Missing residue position could not be inferred: ', count_missing_res)
