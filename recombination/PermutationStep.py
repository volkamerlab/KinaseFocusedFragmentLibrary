
class PermutationStep:

    def __init__(self, fragment, binding_site, depth, subpockets):

        self.fragment = fragment
        self.binding_site = binding_site  # dummy atom
        self.depth = depth  # number of combined fragments
        self.subpockets = subpockets  # list of targeted subpockets
