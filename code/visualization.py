path = '../../data/KLIFS_download/'
radius = 2
from biopandas.mol2 import PandasMol2
from pymol import cmd
# from pymol.cgo import *


def visualizeSubpocketCenters(subpockets, pocket, folder):

    f = open(path+folder+'/subpockets.cgo', 'w')
    # write header
    f.write('# Subpockets of kinase binding pocket\n')
    f.write('from pymol.cgo import *\n')
    f.write('from pymol import cmd\n')
    f.write('obj = [\n')

    # write subpocket centers
    for subpocket in subpockets:
        # transparency
        f.write('\tALPHA, 0.75,\n')
        # color
        f.write('\tCOLOR, '+subpocket.color+',\n')
        # draw sphere
        center = subpocket.center
        f.write('\tSPHERE, %.2f, %.2f, %.2f, %.2f,\n\n' % (center[0], center[1], center[2], radius))

    f.write(']\n')
    f.write('cmd.load_cgo(obj,\'subpockets\')\n')
    f.close()

    return None


def visualizePocket(folder):

    pmol = PandasMol2().read_mol2(path+folder+'/pocket.mol2',
                                  columns={0: ('atom_id', int), 1: ('atom_name', str), 2: ('x', float), 3: ('y', float), 4: ('z', float),
                                           5: ('atom_type', str), 6: ('res_id', int), 7: ('res_name', str), 8: ('charge', float),
                                           9: ('secondary structure', str)})
    cmd.remove('solvent')
    residues = get_bs_res(pmol)

    # color binding site regions
    #for res in residues:
     #   cmd.color()

    return None

# EVA'S FUNCTIONS

def get_bs_res(pmol):
    res_id_list = []
    for Id in pmol.df.res_name:
        if Id[3:] in res_id_list:
            continue
        else:
            res_id_list.append(Id[3:])
    return res_id_list


def klifs_to_pdb(pmol, l):
    pdb_anchors = []
    for elem in l:
        tmp = []
        for id in elem:
            tmp.append(pmol.df[pmol.df.res_id == id]['res_name'].iloc[0][3:])
        pdb_anchors.append(tmp)

    return pdb_anchors
