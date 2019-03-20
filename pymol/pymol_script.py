from biopandas.mol2 import PandasMol2
import pandas as pd
import pymol
# from pymol.cgo import *
from pymol import cmd

import sys
sys.path.append("/home/paula/Masterarbeit/KinaseFocusedFragmentLibrary/code/")
from preprocessing import addMissingResidues, getFolderName, fixResidueIDs
from pocketIdentification import getRegion

path = '/home/paula/Masterarbeit/data/KLIFS_download/'
info = addMissingResidues(pd.read_csv(path + 'overview.csv'))

# Launch PyMol
pymol.finish_launching()

# EVA'S FUNCTIONS

def get_bs_res(pmol):
    res_id_list = []
    klifs_id_list = []
    for i, name in enumerate(pmol.res_name):
        if name[3:] in res_id_list:
            continue
        else:
            res_id_list.append(name[3:])
            klifs_id_list.append(pmol.res_id[i])
    return dict(zip(klifs_id_list, res_id_list))


def klifs_to_pdb(pmol, l):
    pdb_anchors = []
    # regions
    for elem in l:
        tmp = []
        # residues in region
        for id in elem:
            tmp.append(pmol.df[pmol.df.res_id == id]['res_name'].iloc[0][3:])
        pdb_anchors.append(tmp)

    return pdb_anchors


# color binding pocket using PyMOL
def visualizePocket(species, kinase, pdb, alt, chain):

    cmd.reinitialize()

    # overview table for this structure
    global info
    info = info[(info.species == species) & (info.kinase == kinase) & (info.pdb == pdb) & (info.alt == alt) & (info.chain == chain)].iloc[0]

    folder = getFolderName(info)

    # mol2 file of this pocket structure
    pocketMol2 = PandasMol2().read_mol2(path+folder+'/pocket.mol2',
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary structure', str)}).df

    # fix residue IDs
    pocketMol2 = fixResidueIDs(pocketMol2, info.missing_residues)
    # corresponding residue numbers within the protein
    residueDict = get_bs_res(pocketMol2)

    # draw protein
    cmd.load(path+folder+'/pocket.mol2')
    cmd.remove('solvent')
    cmd.show_as('cartoon')
    cmd.color('gray70')
    cmd.set('line_width', 3)
    cmd.set('stick_radius', 0.15)

    # color binding site regions
    for res in residueDict:
        if getRegion(res) == 'hinge':
            cmd.color('deeppurple', 'resi '+residueDict[res])
        elif getRegion(res) == 'DFG':
            cmd.color('tv_blue', 'resi '+residueDict[res])
            if res != 80:
                cmd.show('lines', 'resi ' + residueDict[res])
        elif getRegion(res) == 'g.l':
            cmd.color('forest', 'resi '+residueDict[res])
        elif getRegion(res) == 'alphaC':
            cmd.color('red', 'resi '+residueDict[res])
        elif res == 17:
            cmd.color('yelloworange', 'resi '+residueDict[res])
            cmd.show('lines', 'resi '+residueDict[res])
        elif res == 45:
            cmd.color('sand', 'resi '+residueDict[res])
            cmd.show('lines', 'resi '+residueDict[res])
        elif getRegion(res) == 'linker':
            cmd.color('cyan', 'resi '+residueDict[res])

    # draw subpockets
    exec(open(path+folder+'/subpockets.cgo').read())
    # cmd.load_cgo(obj, 'subpockets')

    # draw ligand
    cmd.load(path+folder+'/ligand.mol2')
    # cmd.select('organic')
    cmd.show_as('sticks', 'org')

    return None


# visualizePocket('Human', 'CDK2', '1h01', 'A', 'A')
visualizePocket('Human', 'EGFR', '1m17', 'A', 'A')

