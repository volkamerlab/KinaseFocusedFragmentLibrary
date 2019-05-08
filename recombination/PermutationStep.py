
class PermutationStep:

    def __init__(self, fragment, dummy, depth, subpockets):

        self.fragment = fragment
        self.dummy = dummy  # dummy atom
        self.depth = depth  # number of combined fragments
        self.subpockets = subpockets  # list of targeted subpockets
