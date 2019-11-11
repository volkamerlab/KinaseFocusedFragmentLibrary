
class Combination:

    """
    Comparable representation of a combination of fragments

    Attributes
    ----------
    frag_ids: frozenset(str)
        Strings representing the fragments that the molecule consists of
    bonds: frozenset(tuple(str))
        Bonds through which the fragments are connected.
        The bonds are stored as tuples of atom IDs.

    Methods
    ----------
    __eq__()
        Two Combination objects are equal if they consist of the same fragments which are connected through the same bonds.

    """

    def __init__(self, frag_ids, bonds=None):
        self.frag_ids = frag_ids
        self.bonds = bonds

    def __eq__(self, other):
        return self.frag_ids == other.frag_ids and self.bonds == other.bonds

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.frag_ids, self.bonds))


class PermutationStep:

    """
    Step of the permutation algorithm representing one single dummy atom of a molecule

    Attributes
    ----------
    compound: Compound object
        the molecule containing the dummy atom
    dummy: str
        atom ID (frag_atom_id) of this atom
    subpocket: str
        Subpocket of the atom adjacent to the dummy atom (current subpocket)
    neighboring_subpocket: str
        Subpocket of the dummy atom
    """

    def __init__(self, mol, dummy, subpocket, neighboring_subpocket):

        self.compound = mol
        self.dummy = dummy
        self.subpocket = subpocket
        self.neighboring_subpocket = neighboring_subpocket


class Compound:

    """
    Represents a combination of fragments including its dummy atoms

    Attributes
    ----------
    frag_ids: list(str)
        Strings representing the fragments that the molecule consists of
    subpockets: list(str)
        Subpockets that the molecule is targeting
    ports: list(Port)
        Port objects representing the dummy atoms of the molecule
    bonds: list(tuple(str))
        Bonds through which the fragments are connected.
        The bonds are stored as tuples of atom IDs.
    """

    def __init__(self, frag_ids, subpockets, ports, bonds):

        self.frag_ids = frag_ids
        self.subpockets = subpockets
        self.ports = ports
        self.bonds = bonds


class Fragment:

    """
    Represents a single fragment from the fragment library

    Attributes
    ----------
    frag_id: str
        ID of the fragment: subpocket_ID, e.g. AP_5
    subpocket: str
        Subpocket that the fragment is targeting
    ports: list(Port)
        Port objects representing the dummy atoms of the fragment
    """

    def __init__(self, frag_id, subpocket, ports):

        self.frag_id = frag_id
        self.subpocket = subpocket  # list of targeted subpockets
        self.ports = ports  # list of Port objects


class Port:

    """
    Represents a single dummy atom

    Attributes
    ----------
    atom_id: str
        frag_atom_id of the dummy atom
    subpocket: str
        Subpocket of the atom adjacent to the dummy atom (subpocket of the fragment containing the dummy)
    neighboring_subpocket: str
        Subpocket of the dummy atom
    bond_type: str
        Type of the bond connecting the dummy to its adjacent atom
    environment: str
        Type of the environment of the current fragment (of the adjacent atom)
    """

    def __init__(self, atom_id, subpocket, neighboring_subpocket, bond_type, environment):

        self.atom_id = atom_id
        self.subpocket = subpocket
        self.neighboring_subpocket = neighboring_subpocket
        self.bond_type = bond_type
        self.environment = environment
