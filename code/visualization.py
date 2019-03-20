path = '../../data/KLIFS_download/'
radius = 2


# create cgo file with subpocket centers as spheres to load in PyMOL
def visualizeSubpocketCenters(subpockets, folder):

    f = open(path+folder+'/subpockets.cgo', 'w')
    # write header
    f.write('# Subpockets of kinase binding pocket\n')
    f.write('from pymol.cgo import *\n')
    f.write('from pymol import cmd\n')
    f.write('obj = [\n')

    # write subpocket centers
    for subpocket in subpockets:
        # transparency
        f.write('\tALPHA, 0.9,\n')
        # color
        f.write('\tCOLOR, '+subpocket.color+',\n')
        # draw sphere
        center = subpocket.center
        f.write('\tSPHERE, %.2f, %.2f, %.2f, %.2f,\n\n' % (center[0], center[1], center[2], radius))

    f.write(']\n')
    f.write('cmd.load_cgo(obj,\'subpockets\')\n')
    f.close()

    return None

