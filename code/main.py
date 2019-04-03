from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2

from pocketIdentification import getSubpocketFromPos, getGeometricCenter, fixSmallFragments, calculateSubpocketCenter
from fragmentation import findBRICSFragments, getFragmentsFromAtomTuples
from classes import Subpocket
from preprocessing import preprocessKLIFSData, getFolderName, getFileName, fixResidueIDs
from visualization import visualizeSubpocketCenters

import os

# ============================= INITIALIZATIONS ===============================================

# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[51], color='0.0, 1.0, 1.0'),  # cyan  # leave out 2? (1m17 will improve)
              Subpocket('AP', residues=[46, 51, 75, 15], color='0.6, 0.1, 0.6'),  # deeppurple
              Subpocket('FP', residues=[72, 51, 4, 81], color='0.2, 0.6, 0.2'),  # forest # substituted 4 for 7 and 72 for 74
              Subpocket('GA', residues=[45, 17, 80], color='1.0, 0.5, 0.0'),  # orange
              # Subpocket('BP', residues=[82, 24, 43], color='0.5, 0.0, 1.0')  # purpleblue
              Subpocket('B1', residues=[81, 28, 43, 38], color='0.0, 0.5, 1.0'),  # marine
              Subpocket('B2', residues=[18, 24, 70, 83], color='0.5, 0.0, 1.0')  # purpleblue # [24, 83, 8, 42] removed 8 because often missing
              ]

# ============================= DATA PREPARATION ============================================

path_to_library = '../FragmentLibrary/'

path_to_data = '../../data/KLIFS_download/'
path_to_KLIFS_download = path_to_data + 'overview.csv'
path_to_KLIFS_export = path_to_data + 'KLIFS_export.csv'

KLIFSData = preprocessKLIFSData(path_to_KLIFS_download, path_to_KLIFS_export)
KLIFSData = KLIFSData[KLIFSData.species == 'Human']
KLIFSData = KLIFSData[KLIFSData.dfg == 'in']
# We are not interested in adenosine phosphates
KLIFSData = KLIFSData[KLIFSData.pdb_id != 'AMP']
KLIFSData = KLIFSData[KLIFSData.pdb_id != 'ADP']
KLIFSData = KLIFSData[KLIFSData.pdb_id != 'ATP']
KLIFSData = KLIFSData[KLIFSData.pdb_id != 'ACP']
KLIFSData = KLIFSData[KLIFSData.pdb_id != 'ANP']
KLIFSData = KLIFSData[(KLIFSData.family == 'RAF') | (KLIFSData.family == 'EGFR') | (KLIFSData.family == 'CDK')]


# clear output files
for subpocket in subpockets:
    folderName = path_to_library+subpocket.name+'/'
    if not os.path.exists(folderName):
        os.mkdir(folderName)
    kinases = set(KLIFSData.kinase)
    for kinase in kinases:
        fileName = folderName + kinase + '.sdf'
        if os.path.isfile(fileName):
            os.remove(fileName)

discardedFragments = []
discardedLigands = []

# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = getFolderName(entry)
    # print(folder, entry.dfg, entry.ac_helix)

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(path_to_data + folder + '/ligand.mol2', removeHs=False)
    pocket = Chem.MolFromMol2File(path_to_data + folder + '/pocket.mol2', removeHs=False)

    try:
        ligandConf = ligand.GetConformer()
    except AttributeError:  # empty molecule
        print('ERROR in ' + folder + ':')
        print('Ligand '+entry.pdb_id+' ('+folder+') could not be loaded. \n')
        continue
    try:
        pocketConf = pocket.GetConformer()
    except AttributeError:
        print('ERROR in ' + folder + ':')
        print('Pocket '+folder+' could not be loaded. \n')
        continue

    # skip multi ligands
    if '.' in Chem.MolToSmiles(ligand):
        print('ERROR in ' + folder + ':')
        print('Ligand consists of multiple molecules. Structure is skipped. \n')
        continue

    lenLigand = ligand.GetNumAtoms()

    # read atom information from binding pocket mol2 file (necessary for residue information)
    # pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')
    pocketMol2 = PandasMol2().read_mol2(path_to_data + folder + '/pocket.mol2',
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary_structure', str)}).df
    # fix residue IDs
    pocketMol2 = fixResidueIDs(pocketMol2, entry.missing_residues)

    # ============================ SUBPOCKET CENTERS =========================================

    skipStructure = False

    # calculate subpocket centers
    for subpocket in subpockets:

        center = calculateSubpocketCenter(subpocket, pocket, pocketMol2, folder)
        # skip structure if no center could be calculated because of missing residues
        if not center:
            skipStructure = True

    # skip this molecule if important residues are missing
    if skipStructure:
        continue

    # visualize subpocket centers using PyMOL
    visualizeSubpocketCenters(subpockets, folder)

    # # get subpocket for each ligand atom
    # for a, atom in enumerate(ligand.GetAtoms()):
    #     # subpocket = getSubpocketFromAtomDistances(a, ligandConf, pocketConf, residues)
    #     subpocket = getSubpocketFromAtom(a, ligandConf, subpockets)
    #     atom.SetProp('subpocket', subpocket)

    # ================================ BRICS FRAGMENTS ==========================================

    skipStructure = False

    # find BRICS fragments and bonds (as atom numbers)
    BRICSFragments, BRICSBonds = findBRICSFragments(ligand)

    # calculate fragment centers and get nearest subpockets
    for BRICSFragment in BRICSFragments:

        center = getGeometricCenter(BRICSFragment.mol.GetAtoms(), BRICSFragment.mol.GetConformer())
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

        subpocket = getSubpocketFromPos(center, subpockets)
        BRICSFragment.subpocket = subpocket

    # discard any fragments where the AP fragment is larger than 20 heavy atoms (e.g. staurosporine)
        if BRICSFragment.mol.GetNumHeavyAtoms() > 22:  # and subpocket == 'AP'
            discardedLigands.append([ligand, entry.pdb])
            skipStructure = True
            break

    if skipStructure:
        continue

    # Deal with small fragments
    fixSmallFragments(BRICSFragments, BRICSBonds)

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
    fragments = getFragmentsFromAtomTuples(bonds, BRICSFragments, ligand)

    # ================================ FRAGMENT LIBRARY ========================================

    # add fragments to their respective pool

    for fragment in fragments:
        # store ligand for this fragment
        fragment.structure = entry.pdb
        # discard large fragments
        if fragment.mol.GetNumHeavyAtoms() > 29:
            discardedFragments.append(fragment)
        else:
            # output_file = path_to_library+fragment.subpocket+'/'+getFileName(entry)+'.sdf'
            output_file = open(path_to_library+fragment.subpocket+'/'+entry.kinase+'.sdf', 'a')
            # print(Chem.MolToMolBlock(fragment.mol), file=open(output_file, 'a'))
            w = Chem.SDWriter(output_file)
            w.write(fragment.mol)

    # ================================ DRAW FRAGMENTS ==========================================

    # convert to 2D molecules for drawing
    for fragment in fragments:
        tmp = AllChem.Compute2DCoords(fragment.mol)

    # atom = fragment.mol.GetAtoms()[0]
    # print(atom.GetProp('atomNumber'), atom.GetProp('neighboringSubpocket'), atom.GetProp('subpocket'), atom.GetProp('priority'),
    #     [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()])

    img = Draw.MolsToGridImage([Chem.RemoveHs(fragment.mol) for fragment in fragments],
                               legends=[fragment.subpocket for fragment in fragments],
                               subImgSize=(400, 400))
    img.save('../fragmented_molecules/'+getFileName(entry)+'.png')


# draw discarded fragments
img = Draw.MolsToGridImage([Chem.RemoveHs(fragment.mol) for fragment in discardedFragments],
                           legends=[fragment.structure+' '+fragment.subpocket+str(fragment.mol.GetNumHeavyAtoms()) for fragment in discardedFragments],
                           subImgSize=(400, 400), molsPerRow=4)
img.save('../discarded_fragments.png')

# draw discarded ligands
for ligand in discardedLigands:
    tmp = AllChem.Compute2DCoords(ligand[0])
img = Draw.MolsToGridImage([Chem.RemoveHs(ligand[0]) for ligand in discardedLigands],
                           legends=[ligand[1] for ligand in discardedLigands],
                           subImgSize=(400, 400))
img.save('../discarded_ligands.png')


# TO DO:
# - store bond information (BRICS rule?)
#
# - deal with missing residues
#
# - deal with multi ligands
