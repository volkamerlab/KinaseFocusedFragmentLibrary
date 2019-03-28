
# LIGAND INFORMATION

ligandConf = ligand.GetConformer()
# for atom in ligand.GetAtoms():
#     # print(atom.GetIdx(), atom.GetAtomicNum(), atom.GetDegree(), [a.GetIdx() for a in atom.GetNeighbors()])
#     idx = atom.GetIdx()
#     pos = ligandConf.GetAtomPosition(idx)
#     print(idx, atom.GetSymbol(), '{0:12.4f}{1:12.4f}{2:12.4f}'.format(pos.x, pos.y, pos.z))


# USING THE ENTIRE COMPLEX AND RDKIT DISTANCE FUNCTION

# complex = Chem.MolFromMol2File('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/complex.mol2', removeHs=False)
# complex = Chem.CombineMols(ligand, pocket)
# complexConf = complex.GetConformer()
# for atom in complex.GetAtoms():
#     idx = atom.GetIdx()
#     pos = complexConf.GetAtomPosition(idx)
#     print(idx, atom.GetSymbol(), '{0:12.4f}{1:12.4f}{2:12.4f}'.format(pos.x, pos.y, pos.z))

# start = time.time()
# distMatrix = pd.DataFrame.from_dict(Chem.Get3DDistanceMatrix(complex))
# #print(distMatrix[5217][2523])
# end = time.time()
# print(end-start)

# complex = Chem.MolFromPDBFile('../../data/KLIFS_download/HUMAN/EGFR/3w2s_altB_chainA/3w2s.pdb')
# res = complex.GetAtoms()[0].GetPDBResidueInfo()
# res_name = res.GetResidueName()
# print(res_name)
