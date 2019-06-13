from pathlib import Path
import pickle
from pickle_loader import pickle_loader


def results_to_file(results, n_results_out):

    """
    Writes the given results (Combination objects) to a pickle and returns the current number of output files

    Parameters
    ----------
    results: set(Combination)
        set of Combination objects representing ligands resulting from the recombinations
    n_results_out: int
        current number of output files containing results

    Returns
    -------
    Int
        Number of output files containing results

    """

    # write results in result set to output file
    n_out = len(results)
    out_path = Path('results/results' + str(n_results_out) + '.pickle')
    with open(out_path, 'wb') as pickle_out:
        print('Write ' + str(n_out) + ' results to file', n_results_out)
        for result in results:
            pickle.dump(result, pickle_out)

    return n_results_out + 1


def add_to_results(combos, results, n_results_out):

    """
    - Checks if the given combos (Combination objects) are already present in the results (including result output files)
    - Adds them to the result set if not

    Parameters
    ----------
    combos: set(Combination)
        - set of Combination objects representing ligands resulting from the recombinations
        - will be added to results
    results: set(Combination)
        set of Combination objects representing ligands resulting from the recombinations
    n_results_out: int
        current number of output files containing results + 1

    """

    # if no results file exists yet, add all new results to result set
    if n_results_out == 0:
        for combo in combos:
            results.add(combo)
        return

    # if results files exist:
    # Check whether new results are already present in these file before adding to result set
    duplicates = set()
    for n in range(n_results_out):

        in_path = Path('results/results' + str(n) + '.pickle')
        with open(in_path, 'rb') as pickle_in:
            results_in = set(pickle_loader(pickle_in))
            for combo in combos:
                if combo in results_in:
                    duplicates.add(combo)

    new_combos = combos - duplicates
    for combo in new_combos:
        results.add(combo)

    return


def process_result(combo, results_temp, results, limit_r, n_results_out, count_results):

    """
    - adds new result to the temporary result set
    - adds results from temp result set to final result set if limit is exceeded
    - writes results to file if limit size is exceeded

    Parameters
    ----------
    combo: Combination
        will be added to results
    results_temp: set(Combination)
        temporary result set
    results: set(Combination)
        final result set
    limit_r: int
        max size of result sets
    n_results_out: int
        current number of output files containing results + 1
    count_results: int
        current number of results
    """

    results_temp.add(combo)
    if len(results_temp) >= limit_r:
        add_to_results(results_temp, results, n_results_out)
        results_temp = set()
    if len(results) >= limit_r:
        count_results += len(results)
        n_results_out = results_to_file(results, n_results_out)
        results = set()

    return results, results_temp, n_results_out