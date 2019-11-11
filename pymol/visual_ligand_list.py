from subpocket_arrangements import subpocket_arrangements
from pymol import cmd

p = '/home/paula/Masterarbeit/'

path = p + 'KinaseFocusedFragmentLibrary/fragmentation/'


def visual_ligand_list(file):

    with open(path+file) as f:

        first = f.readline().strip('\n')
        subpocket_arrangements(first)

        for struct in f.readlines():

            struct = struct.strip('\n')

            cmd.load(p + 'data/KLIFS_download/' + struct + '/ligand.mol2', discrete=1)
            exec(open(p + 'data/KLIFS_download/' + struct + '/subpockets.cgo').read())

    return None

# How to run in PyMol:
# run invalid_subpocket_connections.py
# run subpocket_arrangements.py
# invalid_subpocket_connections('fp-b2')
