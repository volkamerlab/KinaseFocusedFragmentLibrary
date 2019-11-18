radius = 2


def visual_subpockets(subpockets, output_path):

    """
    Create cgo file with subpocket centers as spheres to load in PyMOL

    Parameters
    ----------
    subpockets: list(Subpocket)
        list of Subpocket objects containing the 3D centers and colors of the subpocket
    output_path: Path

    """

    with open(output_path / 'subpockets.cgo', 'w') as f:

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
