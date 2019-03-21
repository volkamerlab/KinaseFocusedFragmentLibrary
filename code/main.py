from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2

from pocketIdentification import getSubpocketFromAtom, getSubpocketFromPos, getGeometricCenter, checkSubpockets, fixSubpockets
from functions import mostCommon, getCaAtom
from fragmentation import FindBRICSFragments, getFragmentsFromAtomTuples
from classes import Fragment, Subpocket
from preprocessing import preprocessKLIFSData, getFolderName, getFileName, fixResidueIDs
from visualization import visualizeSubpocketCenters


# ============================= INITIALIZATIONS ===============================================

# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[50], color='0.0, 1.0, 1.0'),  # cyan  # leave out 2? (1m17 will improve)
              Subpocket('AP', residues=[46, 50, 75, 15], color='0.6, 0.1, 0.6'),  # deeppurple
              Subpocket('FP', residues=[74, 51, 7, 81], color='0.2, 0.6, 0.2'),  # forest # 7 instead of 4  -> but 7 is often missing
              Subpocket('GA', residues=[45, 17, 81], color='1.0, 0.5, 0.0'),  # orange
              Subpocket('BP', residues=[82, 24, 43], color='0.5, 0.0, 1.0')  # purpleblue
              # Subpocket('BP1', residues=[81, 29, 43, 38], color='tv_blue'),
              # Subpocket('BP2', residues=[24, 82, 8], color='forest')
              ]

# ============================= DATA PREPROCESSING ============================================

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
KLIFSData = KLIFSData[KLIFSData.kinase == 'EGFR']


# iterate over molecules
for index, entry in KLIFSData.iterrows():

    # ================================== READ DATA ============================================

    folder = getFolderName(entry)
    print(folder, entry.dfg, entry.ac_helix)

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(path_to_data + folder + '/ligand.mol2', removeHs=False)
    pocket = Chem.MolFromMol2File(path_to_data + folder + '/pocket.mol2', removeHs=False)

    try:
        ligandConf = ligand.GetConformer()
    except AttributeError:  # empty molecule
        print('Ligand '+entry.pdb_id+' ('+folder+') could not be loaded.')
        continue
    try:
        pocketConf = pocket.GetConformer()
    except AttributeError:
        print('Pocket '+folder+' could not be loaded.')
        continue

    # skip multi ligands
    if '.' in Chem.MolToSmiles(ligand):
        print('Ligand consists of multiple molecules. Structure is skipped.')
        continue

    lenLigand = ligand.GetNumAtoms()

    # read atom information from binding pocket mol2 file (necessary for residue information)
    # pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')
    pocketMol2 = PandasMol2().read_mol2(path_to_data + folder + '/pocket.mol2',
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary structure', str)}).df
    # fix residue IDs
    pocketMol2 = fixResidueIDs(pocketMol2, entry.missing_residues)

    # ============================ SUBPOCKET CENTERS =====================================

    skipStructure = False
    # fixSubpockets = False

    # calculate subpocket centers
    for subpocket in subpockets:
        CaAtoms = [getCaAtom(res, pocketMol2, pocket) for res in subpocket.residues]
        # if residue is missing, skip structure
        if None in CaAtoms:
            print('ERROR: Important residue is missing in structure. Structure is skipped.')
            skipStructure = True
            break
            # IDEA: Do not skip structure and use continue instead
            # -> The subpocket center of the previous structure will be used.

        # overwrite subpocket center for current structure
        center = getGeometricCenter(CaAtoms, pocketConf)
        subpocket.center = center
        # print(subpocket.name, subpocket.center)

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

    # ================================ BRICS FRAGMENTS ======================================

    # find BRICS fragments and bonds (as atom numbers)
    BRICSFragments, BRICSBonds = FindBRICSFragments(ligand)

    # calculate fragment centers and get according subpockets
    for BRICSFragment in BRICSFragments:
        center = getGeometricCenter(BRICSFragment.mol.GetAtoms(), BRICSFragment.mol.GetConformer())
        BRICSFragment.center = center
        # calculate distances from subpockets to fragments
        subpocket = getSubpocketFromPos(center, subpockets)
        BRICSFragment.subpocket = subpocket
    # --> Do we still need atom subpockets for other functions (fragmentation?)

    # ================================== FRAGMENTATION ==========================================

    skipStructure = False

    # list to store the bonds where we will cleave (as atom tuples)
    bonds = []
    # BRICS fragments as Fragment objects
    # BRICSFragments = [Fragment(atomNumbers=BRICSFragmentsAtoms[f]) for f in range(len(BRICSFragmentsAtoms))]

    # iterate over BRICS bonds
    for beginAtom, endAtom in BRICSBonds:

        # find corresponding fragments
        firstFragment = [fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers][0]
        secondFragment = [fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers][0]

        # # add subpocket to fragment objects (if not yet defined for this fragment)
        # if firstFragment.subpocket is None:
        #     firstFragment.subpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in firstFragment.atomNumbers])
        #
        # if secondFragment.subpocket is None:
        #     secondFragment.subpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in secondFragment.atomNumbers])

        # # if subpockets of the 2 fragments differ
        # if firstFragment.subpocket != secondFragment.subpocket:
        #
        #     # join very small fragments with larger fragment
        #     if firstFragment.mol.GetNumHeavyAtoms() <= 2:
        #         fixSubpockets(firstFragment, secondFragment)
        #     if secondFragment.mol.GetNumHeavyAtoms() <= 2:
        #         fixSubpockets(secondFragment, firstFragment)

        # check validity of subpockets
        if not checkSubpockets(firstFragment.subpocket, secondFragment.subpocket):

            print("ERROR: Subpockets "+firstFragment.subpocket+" and "+secondFragment.subpocket+" can not be connected."
                                                                                                "Structure is skipped.")
            # skip this molecule if subpocket definition is not valid
            skipStructure = True
            break

        # check if subpockets differ
        if firstFragment.subpocket != secondFragment.subpocket:
            # store this bond as a bond where we will cleave
            bonds.append((beginAtom, endAtom))

    # # iterate over BRICS bonds second time
    # for beginAtom, endAtom in BRICSBonds:
    #
    #     # find corresponding fragments
    #     firstFragment = [fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers][0]
    #     secondFragment = [fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers][0]
    #     # check again if subpockets differ because they may have changed
    #     if firstFragment.subpocket != secondFragment.subpocket:
    #         # store this bond as a bond where we will cleave
    #         bonds.append((beginAtom, endAtom))

    # skip this molecule if subpocket definition is not valid
    if skipStructure:
        continue

    # check if subpockets appear several times and not next to each other

    # # if the ligand is just one BRICS fragment (then the for loop above is not executed because there are no bonds)
    # if len(BRICSFragments) == 1:
    #     subpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in BRICSFragments[0].atomNumbers])
    #     BRICSFragments[0].subpocket = subpocket

    # actual fragmentation
    fragments = getFragmentsFromAtomTuples(bonds, BRICSFragments, ligand)

    # ================================ FRAGMENT LIBRARY ========================================

    # add fragments to their respective pool

    for fragment in fragments:
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


# TO DO:
# - store bond information (BRICS rule?)
#
# - deal with missing residues
#
# - deal with multi ligands
