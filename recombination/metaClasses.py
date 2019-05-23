

class Combination:

    def __init__(self, frag_ids, bonds=None):
        self.frag_ids = frag_ids  # frozenset of fragment IDs
        self.bonds = bonds  # frozenset of tuples of atom indices

    def __eq__(self, other):
        return self.frag_ids == other.frag_ids and self.bonds == other.bonds

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.frag_ids, self.bonds))


class PermutationStep:

    def __init__(self, mol, dummy, subpocket, neighboring_subpocket):

        self.compound = mol  # Compound object
        self.dummy = dummy  # frag atom ID
        self.subpocket = subpocket
        self.neighboring_subpocket = neighboring_subpocket


class Compound:

    def __init__(self, frag_ids, subpockets, ports, bonds):

        # put fragments an subpockets in a dict?
        self.frag_ids = frag_ids  # list of fragment ids
        self.subpockets = subpockets  # list of targeted subpockets
        self.ports = ports  # list of Port objects
        self.bonds = bonds


class Fragment:

    def __init__(self, frag_id, subpocket, ports):

        self.frag_id = frag_id
        self.subpocket = subpocket  # list of targeted subpockets
        self.ports = ports  # list of Port objects


class Port:

    def __init__(self, atom_id, subpocket, neighboring_subpocket):

        self.atom_id = atom_id  # atom ID of this dummy atom
        self.subpocket = subpocket  # subpocket of the current fragment containing the port
        self.neighboring_subpocket = neighboring_subpocket  # neighboring subpocket / subpocket of the dummy atom
