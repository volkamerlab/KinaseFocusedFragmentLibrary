

# load atom information block from a mol2 file into a list of lists of strings
def loadAtomInfoFromMol2(file):

    atomInfo = []
    with open(file) as f:
        line = f.readline()
        while not line.startswith('@<TRIPOS>ATOM'):
            line = f.readline()
        line = f.readline()
        while not line.startswith('@<TRIPOS>'):
            atomInfo.append(line.split())
            line = f.readline()

    return atomInfo


# remove duplicates from array/list x
def removeDuplicates(x):
    return list(dict.fromkeys(x))
