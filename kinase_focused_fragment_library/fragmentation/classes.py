
class Fragment:

    """
    Represents a molecular fragment

    Attributes
    ----------
    mol: Mol
        RDKit molecule representation of the fragment
    atomNumbers: list(int)
        atom numbers of the fragment atoms within the fragmented molecule
    subpocket: Subpocket
        Subpocket in which the fragment lies
    structure: str
        PDB code + alt + chain of corresponding structure
    environment: list[str]
        List of BRICS environment types of the fragment (based on RDKit definitions)
    """

    def __init__(self, mol=None, atomNumbers=None, subpocket=None, structure=None, center=None, environment=None):

        self.mol = mol
        self.atomNumbers = atomNumbers
        self.center = center
        self.subpocket = subpocket
        self.structure = structure
        if environment is None:
            self.environment = []
        elif isinstance(environment, str):
            self.environment = [environment]
        elif isinstance(environment, list):
            self.environment = environment
        else:
            raise ValueError(f"Input type {type(environment)} of `environment` can only be str or list[str]")


class Subpocket:

    """
    Represents a subpocket of the kinase binding site

    Attributes
    ----------
    name: str
    residues: list(int)
        Residues (as KLIFS IDs) used to calculate the subpocket center
    center: list(int)
        3D point representing the geometric center of the subpocket
    color: str
        Color of the subpocket center in PyMol, e.g. '0.0, 1.0, 1.0' (cyan)
    """

    def __init__(self, name, residues=None, center=None, color=None):

        self.name = name
        self.residues = residues
        self.center = center
        self.color = color

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name
