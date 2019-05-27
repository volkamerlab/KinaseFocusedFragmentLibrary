from pathlib import Path
import pickle
from pickle_loader import pickle_loader


def results_to_file(results, n_results_out):

    # write results in result set to output file
    n_out = len(results)
    out_path = Path('results/results' + str(n_results_out) + '.pickle')
    pickle_out = out_path.open('wb')
    print('Write ' + str(n_out) + ' results to file', n_results_out)
    for result in results:
        pickle.dump(result, pickle_out)
    pickle_out.close()

    return n_results_out + 1


def add_to_results(combos, results, n_results_out):

    # if no results file exists yet, add all new results to result set
    if n_results_out == 0:
        for combo in combos:
            results.add(combo)
        return None

    # if results files exist:
    # Check whether new results are already present in these file before adding to result set
    print('Read results from files.')
    duplicates = set()
    for n in range(n_results_out):

        in_path = Path('results/results' + str(n) + '.pickle')
        pickle_in = in_path.open('rb')
        results_in = set(pickle_loader(pickle_in))
        for combo in combos:
            if combo in results_in:
                duplicates.add(combo)
        pickle_in.close()
    new_combos = combos - duplicates
    for combo in new_combos:
        results.add(combo)

    return None


def frags_in_queue_to_file(frags_in_queue, n_frags_out):

    # write objects in frags_in_queue to tmp file
    n_out = len(frags_in_queue)
    out_path = Path('tmp/tmp' + str(n_frags_out) + '.pickle')
    pickle_out = out_path.open('wb')
    print('Write ' + str(n_out) + ' combinations to file', n_frags_out)
    for combo in frags_in_queue:
        pickle.dump(combo, pickle_out)
    pickle_out.close()
    return n_frags_out + 1


def add_to_queue(combos, queue, n_frags_out):

    return False