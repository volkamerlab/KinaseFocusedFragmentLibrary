from biopandas.mol2 import PandasMol2
import pandas as pd
import pymol
# from pymol.cgo import *
from pymol import cmd
from pathlib import Path

import sys
sys.path.append("../fragmentation")
from preprocessing import add_missing_residues, get_folder_name, fix_residue_numbers
from pocketIdentification import get_region

path = '../../data/KLIFS_download/'
info = add_missing_residues(pd.read_csv(path + 'overview.csv'))

# Launch PyMol
pymol.finish_launching()


# function to get dict with 
# keys = KLIFS binding pocket residue numbers 
# values = residue numbers w.r.t. protein
def get_res_dict(pmol):
    resIDs = []
    klifsIDs = []
    for i, name in enumerate(pmol.res_name):
        if name[3:] in resIDs:
            continue
        else:
            resIDs.append(name[3:])
            klifsIDs.append(pmol.res_id[i])
    return dict(zip(klifsIDs, resIDs))


def get_info_from_folder(folder):

    folders = folder.split('/')

    species = folders[0].title()
    kinase = folders[1]
    struct = folders[2]

    if 'alt' in folders[2]:
        struct = struct.split('_')
        pdb = struct[0]
        alt = struct[1][-1]
        chain = struct[2][-1]
    else:
        struct = struct.split('_')
        pdb = struct[0]
        alt = ' '
        chain = struct[1][-1]

    return species, kinase, pdb, alt, chain


# color binding pocket using PyMOL
def subpocket_arrangements(folder):  # (species, kinase, pdb, alt, chain):

    cmd.reinitialize()

    species, kinase, pdb, alt, chain = get_info_from_folder(folder)
    # print(species, kinase, pdb, alt, chain)

    # overview table for this structure
    global info
    info = info[(info.species == species) & (info.kinase == kinase) & (info.pdb == pdb) & (info.alt == alt) & (info.chain == chain)].iloc[0]
    
    # folder = get_folder_name(info)

    # mol2 file of this pocket structure
    pocketMol2 = PandasMol2().read_mol2(path+folder+'/pocket.mol2',
                                        columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                                 5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                                 9: ('secondary structure', str)}).df

    # fix residue IDs
    pocketMol2 = fix_residue_numbers(pocketMol2, info.missing_residues)
    # corresponding residue numbers within the protein
    residueDict = get_res_dict(pocketMol2)

    cmd.bg_color(color='white')
    cmd.set('ray_opaque_background', 'on')
    # draw protein
    cmd.load(path+folder+'/pocket.mol2')
    cmd.remove('solvent')
    cmd.show_as('cartoon')
    cmd.color('yellow')
    cmd.set('line_width', 3)
    cmd.set('stick_radius', 0.15)
    cmd.set('cartoon_transparency', 0.7)

    # color binding site regions
    for res in residueDict:
        if get_region(res) == 'hinge':
            cmd.color('deeppurple', 'resi '+residueDict[res])
        elif get_region(res) == 'DFG':
            cmd.color('tv_blue', 'resi '+residueDict[res])
        elif get_region(res) == 'g.l':
            cmd.color('forest', 'resi '+residueDict[res])
        elif get_region(res).startswith('alpha'):
            cmd.color('red', 'resi '+residueDict[res])
        elif res == 17:
            cmd.color('yelloworange', 'resi '+residueDict[res])
        elif res == 45:
            cmd.color('sand', 'resi '+residueDict[res])
        elif get_region(res) == 'linker':
            cmd.color('cyan', 'resi '+residueDict[res])
        elif get_region(res) == 'b.l':
            cmd.color('forest', 'resi '+residueDict[res])
        elif get_region(res) == 'c.l':
            cmd.color('orange', 'resi '+residueDict[res])

    # draw subpockets
    exec(open(path+folder+'/subpockets.cgo').read())
    # cmd.load_cgo(obj, 'subpockets')

    # draw ligand
    cmd.load(path+folder+'/ligand.mol2')
    # cmd.select('organic')
    cmd.show_as('sticks', 'org')

    return None


# visual_pocket('Human', 'CDK2', '1h01', 'A', 'A')
# visual_pocket('Human', 'EGFR', '1m17', 'A', 'A')
# visual_pocket('Human', 'EGFR', '5em7', 'B', 'A')
# visual_pocket('Human', 'EGFR', '1xkk', ' ', 'A')

