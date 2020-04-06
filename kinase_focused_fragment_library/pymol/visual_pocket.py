from biopandas.mol2 import PandasMol2
import pandas as pd
import pymol
# from pymol.cgo import *
from pymol import cmd
from pathlib import Path

import sys
sys.path.append("../preprocessing")
sys.path.append("../fragmentation")
from preprocessing import add_missing_residues, get_folder_name, fix_residue_numbers
from pocketIdentification import get_region

path = Path('..') / '..'/ '..' / 'KinaseFocusedFragmentLibraryData' / 'KLIFS_download/'
info = add_missing_residues(pd.read_csv(path / 'overview.csv'))

# Launch PyMol
pymol.finish_launching()
cmd.reinitialize()


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


def visual_pocket(species, kinase, pdb, alt, chain, state=1):
    """
    Color binding pocket using PyMOL.
    """

    # overview table for this structure
    global info
    info_selected = info[(info.species == species) & (info.kinase == kinase) & (info.pdb == pdb) & (info.alt == alt) & (info.chain == chain)]

    try:
        info_selected = info_selected.iloc[0]
    except IndexError:
        print(f'No or more than one entry/entries was/were selected, only one allowed.')
        print(info_selected)
    
    folder = get_folder_name(info_selected)

    # mol2 file of this pocket structure
    pocket_mol2 = PandasMol2().read_mol2(
        str(path / folder / 'pocket.mol2'),
        columns={
            0: ('atom_id', int),
            1: ('atom_name', str),
            2: ('x', float),
            3: ('y', float),
            4: ('z', float),
            5: ('atom_type', str),
            6: ('res_id', int),
            7: ('res_name', str),
            8: ('charge', float),
            9: ('secondary structure', str)
        }
    ).df

    # fix residue IDs
    pocket_mol2 = fix_residue_numbers(pocket_mol2, info_selected.missing_residues)

    # corresponding residue numbers within the protein
    residue_dict = get_res_dict(pocket_mol2)

    # draw protein
    cmd.load(str(path / folder / 'pocket.mol2'), 'pocket', state)
    cmd.select(name=f'pocket{state}', selection='pocket', state=state)
    cmd.remove(f'pocket{state} and solvent')

    cmd.show_as('cartoon', f'pocket{state}')
    cmd.color('gray70', f'pocket{state}')
    cmd.set('line_width', 3, f'pocket{state}')
    cmd.set('stick_radius', 0.15, f'pocket{state}')

    # color binding site regions
    for res in residue_dict:
        if get_region(res) == 'hinge':
            cmd.color('deeppurple', f'pocket{state} and resi {residue_dict[res]}')
        elif get_region(res) == 'DFG':
            cmd.color('tv_blue', f'pocket{state} and resi {residue_dict[res]}')
            if res != 80:
                cmd.show('lines', f'pocket{state} and resi {residue_dict[res]}')
        elif get_region(res) == 'g.l':
            cmd.color('forest', f'pocket{state} and resi {residue_dict[res]}')
        elif get_region(res) == 'alphaC':
            cmd.color('red', f'pocket{state} and resi {residue_dict[res]}')
        elif res == 17:
            cmd.color('yelloworange', f'pocket{state} and resi {residue_dict[res]}')
            cmd.show('sticks', f'pocket{state} and resi {residue_dict[res]}')
        elif res == 45:
            cmd.color('sand', f'pocket{state} and resi {residue_dict[res]}')
            cmd.show('lines', f'pocket{state} and resi {residue_dict[res]}')
        elif get_region(res) == 'linker':
            cmd.color('cyan', f'pocket{state} and resi {residue_dict[res]}')

    cmd.show_as('sticks', f'pocket{state} and resi {residue_dict[24]}')
    cmd.delete(f'pocket{state}')

    # draw ligand
    cmd.load(str(path / folder / 'ligand.mol2'), 'ligand', state)
    cmd.select(name=f'ligand{state}', selection='ligand', state=state)
    cmd.show_as('sticks', f'ligand{state} and org')
    cmd.delete(f'ligand{state}')

    # draw subpockets
    try:
        exec(open(path / folder / 'subpockets.cgo').read())
    except FileNotFoundError:
        print(
            f'Subpocket file not found for state {state}: {species}|{kinase}|{pdb}|{alt}|{chain}')


def visual_pockets(n=20):

    # overview table for all structures
    global info

    # Reset rows for PyMol states
    info.reset_index(inplace=True)

    for index, row in info.head(n).iterrows():

        visual_pocket(row.species, row.kinase, row.pdb, row.alt, row.chain, state=index+1)


# visual_pocket('Human', 'CDK2', '1h01', 'A', 'A')
# visual_pocket('Human', 'EGFR', '1m17', 'A', 'A')
# visual_pocket('Human', 'EGFR', '5em7', 'B', 'A')
# visual_pocket('Human', 'EGFR', '1xkk', ' ', 'A')

