
class PermutationStep:

    """
    Step of the permutation algorithm representing one single dummy atom of a molecule

    Attributes
    ----------
    fragment: RDKit Mol object

    dummy: int
        A dummy atom of the fragment (as atom index)

    depth: int
        Number of original fragments that this fragment consists of
    subpockets: list(string)
        List of subpockets that the fragment is targeting

    """

    def __init__(self, fragment, dummy, depth, subpockets):

        self.fragment = fragment
        self.dummy = dummy
        self.depth = depth
        self.subpockets = subpockets
