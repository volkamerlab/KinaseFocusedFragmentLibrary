
class Fragment:

    def __init__(self, mol=None, atomNumbers=None, subpocket=None, smiles=None, ligand=None):

        self.mol = mol  # rdkit molecule representation of the fragment
        self.atomNumbers = atomNumbers  # atom numbers of the fragment atoms within the fragmented molecule
        self.subpocket = subpocket
        self.smiles = smiles
        self.ligand = ligand
        # currently, ligand is a molecule -- instead store the ligand name or so?


class Subpocket:

    def __init__(self, name, residues=None, center=None, color=None):

        self.name = name
        self.residues = residues    # residues used to calculate the subpocket center
        self.center = center
        self.color = color

