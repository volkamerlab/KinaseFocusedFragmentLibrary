import threading

from pickle_loader import pickle_loader
from analyze_results import analyze_result


class LigandAnalysisThread(threading.Thread):

    def __init__(self, files, data, results=None):

        threading.Thread.__init__(self)
        self.files = files
        self.data = data
        self.results = results

    def run(self):

        self.results = []

        for file in self.files:

            print(str(file))
            with open(file, 'rb') as pickle_in:

                for meta in pickle_loader(pickle_in):

                    result = analyze_result(meta, self.data)

                    if result is not None:
                        self.results.append(result)
