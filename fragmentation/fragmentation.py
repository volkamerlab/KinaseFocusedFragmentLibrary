from rdkit import Chem
from rdkit.Chem import BRICS

from classes import Fragment


def find_brics_fragments(mol):

    """
    Carries out the BRICS fragmentation algorithm and returns the BRICS fragments and bonds

    Parameters
    ----------
    mol: RDKit Mol object
        molecule to be fragmented

    Returns
    -------
    fragments: list(Fragment)
        list of the resulting BRICS fragments as Fragment objects
    brics_bonds: list(tuple(tuple(int), tuple(str)))
        list of tuples, where each tuple represents a BRICS bond between two atoms in the ligand
        - first tuple: atom indices
        - second tuple: BRICS environment types
    """

    brics_bonds = list(BRICS.FindBRICSBonds(mol))  # convert generator to list for later purposes
    atom_tuples = [bond[0] for bond in brics_bonds]

    # if mol was not fragmented:
    if len(atom_tuples) == 0:
        fragments = [Fragment(atomNumbers=range(mol.GetNumAtoms()), mol=mol, environment='na')]
        return fragments, atom_tuples
    # else:
    bonds = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atom_tuples]
    broken_mol = Chem.FragmentOnBonds(mol, bonds)  # addDummies=False) # dummy atoms are needed since storing the environment types

    fragment_atoms = Chem.GetMolFrags(broken_mol)
    fragment_mols = Chem.GetMolFrags(broken_mol, asMols=True)

    fragments = [Fragment(atomNumbers=n, mol=m) for (n, m) in zip(fragment_atoms, fragment_mols)]

    return fragments, brics_bonds


def fragmentation(ligand, atom_tuples, brics_fragments):

    """
    - Carries out a fragmentation of the given molecule at the given bonds (which should be a subsection of its BRICS bonds)
    - Assigns each resulting fragment to a subpocket based on the subpocket of the BRICS fragments that it consists of
    - Assigns each atom of the new fragments to a subpocket -> Neighboring subpocket is stored at the dummy atoms

    Parameters
    ----------
    ligand: RDKit Mol object
        molecule to be fragmented
    atom_tuples: list(tuple(int))
            list of atom index tuples, where each tuple represents a bond between two atoms in the ligand
    brics_fragments: list(Fragment)
        list of BRICS fragments of the ligand as Fragment objects

    Returns
    -------
    fragments: list(Fragment)
        list of the resulting fragments as Fragment objects

    """

    # get rdkit bonds (NOT BRICS bonds but custom bonds already)
    bonds = [ligand.GetBondBetweenAtoms(x, y).GetIdx() for x, y in atom_tuples]
    if len(bonds) > 0:
        # fragment ligand at bonds and keep dummy atoms
        fragmented_ligand = Chem.FragmentOnBonds(ligand, bonds)
    else:
        fragmented_ligand = ligand
    # get smiles of fragments
    fragment_smiles = Chem.MolToSmiles(fragmented_ligand).split('.')
    # get rdkit molecules of fragments
    fragment_mols = Chem.GetMolFrags(fragmented_ligand, asMols=True)
    # get atom numbers (w.r.t. ligand) of fragments
    fragment_atoms = Chem.GetMolFrags(fragmented_ligand)

    fragments = []

    # iterate over new fragments
    for (atomNumbers, mol, smile) in zip(fragment_atoms, fragment_mols, fragment_smiles):

        # get subpocket corresponding to fragment (Is there a better way?)
        subpocket = next(brics_fragment.subpocket for brics_fragment in brics_fragments if atomNumbers[0] in brics_fragment.atomNumbers)
        # create Fragment object
        fragment = Fragment(mol=mol, smiles=smile, atomNumbers=atomNumbers, subpocket=subpocket)

        # set atom properties for the created fragment
        for atom, atomNumber in zip(fragment.mol.GetAtoms(), fragment.atomNumbers):

            # get environment type of the brics fragment that the current atom belongs to
            env_type = next(brics_fragment.environment for brics_fragment in brics_fragments if atomNumber in brics_fragment.atomNumbers)

            # if atom is not a dummy atom
            if atom.GetSymbol() != '*':
                # set atom number within the entire molecule as property of the fragment atom
                # IS THIS ALWAYS TRUE? (Does order of atoms always stay the same after fragmentation?)
                atom.SetIntProp('atomNumber', atomNumber)
                atom.SetProp('subpocket', subpocket)
                atom.SetProp('environment', env_type)

            # if atom = dummy atom
            else:
                # neighbor = atom next to a dummy (Can several neighbors exist?)
                neighbor = atom.GetNeighbors()[0]

                # -> This works only because dummy atoms are always last in the iteration
                neighbor_atom = neighbor.GetIntProp('atomNumber')
                # get and set atom number w.r.t ligand of the dummy atom
                bond_atoms = next(atomTuple for atomTuple in atom_tuples if neighbor_atom in atomTuple)
                dummy_atom = next(atomNumber for atomNumber in bond_atoms if atomNumber != neighbor_atom)
                atom.SetIntProp('atomNumber', dummy_atom)
                # get and set subpocket of the dummy atom
                neighboring_subpocket = next(BRICSFragment.subpocket for BRICSFragment in brics_fragments
                                             if dummy_atom in BRICSFragment.atomNumbers)
                atom.SetProp('subpocket', neighboring_subpocket)
                atom.SetProp('environment', env_type)

        fragments.append(fragment)

        # solve linker problem here?
        # 1. get all brics fragments for this fragment
        # 2. identify ambiguous brics fragment(s)
        # 3. remove this substructure from current fragment
        # 4. store resulting fragment additionally, set 'stripped' property?
        # What if resulting fragment is too small? -> check for that; do not strip
        # Will atom properties be kept after editing molecule?

    return fragments
