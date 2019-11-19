from classes_meta import PermutationStep


def add_to_queue(queue, compound):

    """
    Adds each Port of a given Compound object to the queue if the subpocket of the Port is not occupied by the compound yet

    Parameters
    ----------
    queue: deque(PermutationStep)
    compound: Compound

    Returns
    -------
    True if one or multiple Ports have been added to the queue,
    False otherwise

    """

    added_to_queue = False
    for port in compound.ports:
        if port.neighboring_subpocket in compound.subpockets:
            continue

        new_ps = PermutationStep(mol=compound, dummy=port.atom_id, subpocket=port.subpocket,
                                 neighboring_subpocket=port.neighboring_subpocket)
        queue.append(new_ps)
        added_to_queue = True

    return added_to_queue
