from rdkit import Chem

from construct_ligand import construct_ligand
from drug_likeliness import is_drug_like


def analyze_result(meta, data):

    ligand = construct_ligand(meta, data)
    # if ligand could not be constructed, skip
    if not ligand:
        return

    count_ligands += 1

    # number of occupied subpockets
    n_sp[len(meta.frag_ids)] += 1
    # occupied subpockets
    for frag_id in meta.frag_ids:
        n_per_sp[frag_id[:2]] += 1

    # store ligand in combinatorial library
    # property_ligand = PropertyMol(ligand)
    # pickle.dump(property_ligand, combinatorial_library_file)
    # ligand_smiles.add(Chem.MolToSmiles(ligand))

    # necessary for Lipinski rule
    # Chem.SanitizeMol(ligand)  # slower
    Chem.GetSymmSSSR(ligand)

    # Lipinski rule
    lipinski, wt, logp, hbd, hba = is_drug_like(ligand)
    # recombined ligand fingerprint
    # fp = [1] if lipinski else [0]
    if lipinski:
        lipinski_ligands += 1
    wt_ligands += wt
    logp_ligands += logp
    hbd_ligands += hbd
    hba_ligands += hba
    fp.extend([wt, logp, hbd, hba])

    # PAINS substructure search
    match = pains.GetFirstMatch(ligand)
    # if pains was found
    if match is not None:
        count_pains += 1
    # if no pains was found AND Lipinski rule fulfilled
    else:
        if lipinski:
            n_filtered_sp[len(meta.frag_ids)] += 1
            n_atoms_filtered[n] = n_atoms_filtered[n] + 1 if n in n_atoms_filtered else 1
            for frag_id in meta.frag_ids:
                n_filtered_per_sp[frag_id[:2]] += 1
            filtered_ligands += 1

    # number of atoms
    n = ligand.GetNumHeavyAtoms()
    n_atoms[n] = n_atoms[n] + 1 if n in n_atoms else 1

    return  # fp
