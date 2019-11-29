

class Result:

    def __init__(self, meta, lipinski, mwt, logp, hbd, hba, n_atoms, original, original_sub, chembl_match, scaffold):

        self.meta = meta
        self.n_subpockets = len(meta.frag_ids)

        # boolean values
        self.lipinski = lipinski
        self.mwt = mwt
        self.logp = logp
        self.hbd = hbd
        self.hba = hba
        self.n_atoms = n_atoms
        self.original = original
        self.original_sub = original_sub
        self.chembl_match = chembl_match
        self.scaffold = scaffold
