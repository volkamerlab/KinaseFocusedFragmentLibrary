
class Fragment:

    def __init__(self, mol=None, atomNumbers=None, subpocket=None, smiles=None, structure=None, center=None):

        self.mol = mol  # rdkit molecule representation of the fragment
        self.atomNumbers = atomNumbers  # atom numbers of the fragment atoms within the fragmented molecule
        self.center = center
        self.subpocket = subpocket
        self.smiles = smiles
        self.structure = structure  # PDB code + alt + chain of corresponding structure


class Subpocket:

    def __init__(self, name, residues=None, center=None, color=None):

        self.name = name
        self.residues = residues    # residues used to calculate the subpocket center
        self.center = center
        self.color = color

