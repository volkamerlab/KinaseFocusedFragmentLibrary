from rdkit import Chem


def queue_to_tmp(queue, limit, tmp_q_path):

    n_out = int(3/4 * limit)
    tmp_q_file = tmp_q_path.open('a')
    print('Write ' + str(n_out) + ' queue objects to file.')
    w = Chem.SDWriter(tmp_q_file)
    for i in range(n_out):
        ps_out = queue.pop()
        mol = ps_out.fragment
        mol.SetIntProp('dummy', ps_out.dummy)
        mol.SetIntProp('depth', ps_out.depth)
        mol.SetProp('subpockets', ' '.join(ps_out.subpockets))
        Chem.CreateAtomStringPropertyList(mol, 'subpocket')
        try:
            w.write(mol)
        except ValueError as e:
            print('ERROR: ', e)
            continue
    w.flush()
    tmp_q_file.flush()
