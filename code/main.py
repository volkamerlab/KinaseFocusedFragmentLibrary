from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from biopandas.mol2 import PandasMol2
import time

from pocketIdentification import getSubpocketFromAtom, checkSubpockets
from functions import mostCommon, getGeometricCenter, getCaAtom
from fragmentation import FindBRICSFragments, getFragmentsFromAtomTuples
from classes import Fragment, Subpocket
from preprocessing import preprocessKLIFSData, getFolderName

# data preprocessing
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
# KLIFSData = KLIFSData[KLIFSData.kinase == 'EGFR']


# define the 6 subpockets
subpockets = [Subpocket('SE', residues=[50]),
              Subpocket('AP', residues=[48, 12, 77, 45]),
              Subpocket('FP', residues=[74, 51, 4, 81]),
              Subpocket('GA', residues=[45, 17]),
              Subpocket('BP1', residues=[81, 29, 43, 38]),
              Subpocket('BP2', residues=[24, 82, 8])]


# iterate over molecules
for index, entry in KLIFSData.iterrows():

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

    residues = pocketMol2.res_id.apply(int)

    start = time.time()

    # calculate subpocket centers
    for subpocket in subpockets:
        # Do this (get ca atoms) only once for each kinase?
        CaAtoms = [getCaAtom(res, pocketMol2, pocket) for res in subpocket.residues]
        # overwrite subpocket center for current structure
        center = getGeometricCenter(CaAtoms, pocketConf)
        subpocket.center = center
        print(subpocket.name, subpocket.center)

    # get subpocket for each ligand atom
    # TO DO: use subpocket centers for this!
    for a, atom in enumerate(ligand.GetAtoms()):

        subpocket = getSubpocketFromAtom(a, ligandConf, pocketConf, residues)
        atom.SetProp('subpocket', subpocket)

    end = time.time()
    print("Subpocket identification:", end - start)
    start = time.time()

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

    # actual fragmentation
    fragments = getFragmentsFromAtomTuples(bonds, BRICSFragments, ligand)

    for fragment in fragments:
        tmp = AllChem.Compute2DCoords(fragment.mol)

    # atom = fragment.mol.GetAtoms()[0]
    # print(atom.GetProp('atomNumber'), atom.GetProp('neighboringSubpocket'), atom.GetProp('subpocket'), atom.GetProp('priority'),
    #     [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()])

    img = Draw.MolsToGridImage([fragment.mol for fragment in fragments],
                               legends=[fragment.subpocket for fragment in fragments],
                               subImgSize=(400, 400))
    img.save('test/'+entry.pdb+'_subpocket-fragmentation.png')

    end = time.time()
    print("Fragmentation:", end - start)

    # add fragments to their respective pool


# TO DO:
# - store bond information (BRICS rule?)
#
# - implement correct subpocket definition (subpocket centers)
