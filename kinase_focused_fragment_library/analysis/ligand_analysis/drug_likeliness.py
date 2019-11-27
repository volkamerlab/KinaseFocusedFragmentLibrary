from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski


# takes about 1s for 2000 mols
def is_drug_like(mol):

    mol_wt = 1 if Descriptors.ExactMolWt(mol) <= 500 else 0

    logp = 1 if Descriptors.MolLogP(mol) <= 5 else 0

    hbd = 1 if Lipinski.NumHDonors(mol) <= 5 else 0

    hba = 1 if Lipinski.NumHAcceptors(mol) <= 10 else 0

    lipinski = 1 if mol_wt + logp + hbd + hba >= 3 else 0

    return lipinski, mol_wt, logp, hbd, hba
