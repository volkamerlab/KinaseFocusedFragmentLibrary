from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2

from pocketIdentification import getSubpocketFromAtom, getGeometricCenter, checkSubpockets
from functions import mostCommon, getCaAtom
from fragmentation import FindBRICSFragments, getFragmentsFromAtomTuples
from classes import Fragment, Subpocket
from preprocessing import preprocessKLIFSData, getFolderName
from visualization import visualizeSubpocketCenters


# ============================= INITIALIZATIONS ===============================================

# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[50, 2], color='0.0, 1.0, 1.0'),  # cyan
              Subpocket('AP', residues=[46, 50, 75, 12], color='0.6, 0.1, 0.6'),  # deeppurple
              Subpocket('FP', residues=[74, 51, 4, 81], color='1.0, 0.5, 1.0'),  # violet
              Subpocket('GA', residues=[45, 17, 81], color='1.0, 0.5, 0.0'),  # orange
              Subpocket('BP', residues=[82, 24, 43], color='0.3, 0.3, 1.0')  # tv_blue
              # Subpocket('BP1', residues=[81, 29, 43, 38], color='tv_blue'),
              # Subpocket('BP2', residues=[24, 82, 8], color='forest')
              ]

# ============================= DATA PREPROCESSING ============================================

path = '../../data/KLIFS_download/'
path_to_KLIFS_download = path+'overview.csv'
path_to_KLIFS_export = path+'KLIFS_export.csv'

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

    skip_molecule = False

    # load ligand and binding pocket to rdkit molecules
    ligand = Chem.MolFromMol2File(path+folder+'/ligand.mol2', removeHs=False)
    pocket = Chem.MolFromMol2File(path+folder+'/pocket.mol2', removeHs=False)

    try:
        # get molecule conformers
        ligandConf = ligand.GetConformer()
    except:
        print('Ligand '+entry.pdb_id+' ('+folder+') could not be loaded.')
        continue
    try:
        pocketConf = pocket.GetConformer()
    except:
        print('Pocket '+folder+' could not be loaded.')
        continue

    lenLigand = ligand.GetNumAtoms()

    # read atom information from binding pocket mol2 file (necessary for residue information)
    # pocketMol2 = loadAtomInfoFromMol2('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altA_chainA/pocket.mol2')
    pocketMol2 = PandasMol2().read_mol2(path+folder+'/pocket.mol2',
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary structure', str)}).df
    # fix residue IDs
    for res in entry.missing_residues:
        pocketMol2['res_id'].mask(pocketMol2['res_id'] >= res, pocketMol2['res_id']+1, inplace=True)

    # ============================ SUBPOCKET IDENTIFICATION =====================================

    # calculate subpocket centers
    for subpocket in subpockets:
        # Do this (get ca atoms) only once for each kinase?
        CaAtoms = [getCaAtom(res, pocketMol2, pocket) for res in subpocket.residues]
        if None in CaAtoms:
            skip_molecule = True
            break

        # overwrite subpocket center for current structure
        center = getGeometricCenter(CaAtoms, pocketConf)
        subpocket.center = center
        # print(subpocket.name, subpocket.center)

    # skip this molecule if important residues are missing
    if skip_molecule:
        continue

    # visualize subpocket centers using PyMOL
    visualizeSubpocketCenters(subpockets, pocket, folder)

    # get subpocket for each ligand atom
    for a, atom in enumerate(ligand.GetAtoms()):
        # subpocket = getSubpocketFromAtomDistances(a, ligandConf, pocketConf, residues)
        subpocket = getSubpocketFromAtom(a, ligandConf, subpockets)
        atom.SetProp('subpocket', subpocket)

    # ================================== FRAGMENTATION ==========================================

    # find BRICS fragments and bonds (as atom numbers)
    BRICSFragmentsAtoms, BRICSBonds = FindBRICSFragments(ligand)

    # list to store the bonds where we will cleave (as atom tuples)
    bonds = []
    # BRICS fragments as Fragment objects
    BRICSFragments = [Fragment(atomNumbers=BRICSFragmentsAtoms[f]) for f in range(len(BRICSFragmentsAtoms))]

    # iterate over BRICS bonds
    for beginAtom, endAtom in BRICSBonds:

        # find corresponding fragments
        firstFragment = [fragment for fragment in BRICSFragments if beginAtom in fragment.atomNumbers][0]
        secondFragment = [fragment for fragment in BRICSFragments if endAtom in fragment.atomNumbers][0]

        # add subpocket to fragment objects (if not yet defined for this fragment)
        if firstFragment.subpocket is None:
            firstSubpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in firstFragment.atomNumbers])
            firstFragment.subpocket = firstSubpocket
        else:
            firstSubpocket = firstFragment.subpocket
        if secondFragment.subpocket is None:
            secondSubpocket = mostCommon([ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in secondFragment.atomNumbers])
            secondFragment.subpocket = secondSubpocket
        else:
            secondSubpocket = secondFragment.subpocket

        # check validity of subpockets
        if not checkSubpockets(firstSubpocket, secondSubpocket):

            print("ERROR: Subpockets "+firstSubpocket+" and "+secondSubpocket+" can not be connected. "
                                                                              "Molecule is skipped.")
            # skip this molecule if subpocket definition is not valid
            skip_molecule = True
            break

        # if subpockets of the 2 fragments differ (and they are valid)
        if firstSubpocket != secondSubpocket:

            # store this bond as a bond where we will cleave
            bonds.append((beginAtom, endAtom))

    # skip this molecule if subpocket definition is not valid
    if skip_molecule:
        continue

    # if the ligand is just one BRICS fragment (then the for loop above is not executed because there are no bonds)
    if len(BRICSFragments) == 1:
        subpocket = [ligand.GetAtomWithIdx(a).GetProp('subpocket') for a in BRICSFragments[0].atomNumbers]
        BRICSFragments[0].subpocket = subpocket

    # actual fragmentation
    fragments = getFragmentsFromAtomTuples(bonds, BRICSFragments, ligand)

    # ================================ DRAW FRAGMENTS ==========================================

    for fragment in fragments:
        tmp = AllChem.Compute2DCoords(fragment.mol)

    # atom = fragment.mol.GetAtoms()[0]
    # print(atom.GetProp('atomNumber'), atom.GetProp('neighboringSubpocket'), atom.GetProp('subpocket'), atom.GetProp('priority'),
    #     [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()])

    try:

        img = Draw.MolsToGridImage([fragment.mol for fragment in fragments],
                                   legends=[fragment.subpocket for fragment in fragments],
                                   subImgSize=(400, 400))
        img.save('../fragmented_molecules/'+entry.pdb+'.png')
    except:
        print('ERROR: Image could not be created.')

    # ================================ FRAGMENT LIBRARY ========================================

    # add fragments to their respective pool


# TO DO:
# - store bond information (BRICS rule?)
#
# - implement correct subpocket definition (subpocket centers)
