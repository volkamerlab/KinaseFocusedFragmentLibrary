
class PermutationStep:

    """
    Step of the permutation algorithm representing one single dummy atom of a molecule

    Attributes
    ----------
    fragment: RDKit Mol object

    dummy: RDKit Atom object
        A dummy atom of the fragment

    depth: int
        Number of original fragments that this fragment consists of
    subpockets: list(string)
        List of subpockets that the fragment is targeting

    """

    def __init__(self, fragment, dummy, depth, subpockets):

        self.fragment = fragment
        self.dummy = dummy  # dummy atom
        self.depth = depth  # number of combined fragments
        self.subpockets = subpockets  # list of targeted subpockets
