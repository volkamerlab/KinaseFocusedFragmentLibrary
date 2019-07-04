from classes_meta import PermutationStep


def add_to_queue(queue, compound):

    added_to_queue = False
    for port in compound.ports:
        if port.neighboring_subpocket in compound.subpockets:
            continue

        new_ps = PermutationStep(mol=compound, dummy=port.atom_id, subpocket=port.subpocket,
                                 neighboring_subpocket=port.neighboring_subpocket)
        queue.append(new_ps)
        added_to_queue = True

    return added_to_queue
