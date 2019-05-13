import itertools
import pickle


def slice_deque(d, start, stop):
    d.rotate(-start)
    sliced = list(itertools.islice(d, 0, stop-start, 1))
    d.rotate(start)
    return sliced


def pickle_loader(pickle_file):
    try:
        while True:
            yield pickle.load(pickle_file)
    except EOFError:
        pass
