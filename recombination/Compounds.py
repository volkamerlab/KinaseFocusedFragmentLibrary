
class Compound:

    def __init__(self, pdbs, subpockets, ports):

        # put fragments an subpockets in a dict?
        self.pdbs = pdbs  # list of fragment ids
        self.subpockets = subpockets  # list of targeted subpockets
        self.ports = ports  # list of Port objects


class Fragment:

    def __init__(self, pdb, subpocket, ports):

        self.pdb = pdb
        self.subpocket = subpocket  # list of targeted subpockets
        self.ports = ports  # list of Port objects


class Port:

    def __init__(self, subpocket, neighboring_subpocket, fragment=None, dummy=None):

        self.subpocket = subpocket  # subpocket of the current fragment containing the port
        self.neighboring_subpocket = neighboring_subpocket  # neighboring subpocket / subpocket of the dummy atom
        self.fragment = fragment  # Fragment to replace the dummy atom of the port
        self.dummy = dummy  # Dummy atom of the fragment connecting to this port
