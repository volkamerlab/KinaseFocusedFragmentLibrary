import pickle


def pickle_loader(pickle_file):

    """
    Load a pickle file with multiple objects

    Parameters
    ----------
    pickle_file: file object
        input binary pickle file

    Returns
    -------
    Generator object with loaded objects

    """

    try:
        while True:
            yield pickle.load(pickle_file)
    except EOFError:
        pass
