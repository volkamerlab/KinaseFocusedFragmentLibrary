
class Fragment:

    def __init__(self, mol=None, atomNumbers=None, subpocket=None, smiles=None):

        self.mol = mol  # rdkit molecule representation of the fragment
        self.atomNumbers = atomNumbers  # atom numbers of the fragment atoms within the fragmented molecule
        self.subpocket = subpocket
        self.smiles = smiles
