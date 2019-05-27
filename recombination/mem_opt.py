from pathlib import Path
import pickle
from pickle_loader import pickle_loader


def results_to_file(results, n_results_out):

    n_out = len(results)
    out_path = Path('results/results' + str(n_results_out) + '.pickle')
    pickle_out = out_path.open('wb')
    print('Write ' + str(n_out) + ' results to file', n_results_out)
    for result in results:
        pickle.dump(result, pickle_out)
    pickle_out.close()

    return n_results_out + 1


def is_in_results(combo, n_results_out):

    for n in range(n_results_out):

        print('Read results from file', n)
        in_path = Path('results/results' + str(n) + '.pickle')
        pickle_in = in_path.open('rb')
        try:
            while True:
                result = pickle.load(pickle_in)
                if result == combo:
                    return True
        except EOFError:
            pass
        pickle_in.close()

    return False


def frags_in_queue_to_file(frags_in_queue, n_frags_out):

    return n_frags_out + 1


def is_in_queue(combo, n_frags_out):

    return False