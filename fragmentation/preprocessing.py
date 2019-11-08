import pandas as pd


# code taken and adapted from
# https://github.com/AndreaVolkamer/KinaseSimilarity/blob/master/Bachelorarbeit/structFP_phchem/structFP_phchem_preprocessing.ipynb

def preprocess_klifs_data(path_to_klifs_download, path_to_klifs_export):

    """
    Given the two CSV files downloaded from KLIFS:
    - merges the two dataframes
    - adds a column listing the missing residue numbers
    - for each PDB code, selects the structure with the best quality score
    - selects relevant columns
    - returns the filtered dataframe

    Parameters
    ----------
    path_to_klifs_download: Path or Str
        path to overview.csv file given by KLIFS (including the file name)
    path_to_klifs_export: Path or Str
        path KLIFS_export.csv file (including the file name)

    Returns
    -------
    df_screened: Pandas DataFrame
        Filtered dataframe including best quality structure for each PDB

    """

    # read overview file
    df = pd.read_csv(path_to_klifs_download)
    # read export file
    df_csv = pd.read_csv(path_to_klifs_export)
    # rename columns
    df_csv.columns = ['kinase', 'family', 'group', 'pdb', 'chain', 'alt', 'species', 'ligand', 'pdb_id',
                      'allosteric_name', 'allosteric_PDB', 'dfg', 'ac_helix']

    # sync kinase names for both data frames
    for ix, row in df_csv.iterrows():
        a = row.alt
        s = row.kinase
        if '(' in s:
            row.kinase = s[s.find("(") + 1:s.find(")")]
        if a == '-':
            row.alt = ' '

    # outer merge, as information from both data frames should be kept
    df_screen = pd.merge(df, df_csv, how='outer', on=['species', 'kinase', 'pdb', 'chain', 'alt', 'allosteric_PDB'])

    # loose irrelevant data
    df_screen = df_screen[['kinase', 'family', 'group', 'species', 'pdb', 'pdb_id', 'alt', 'chain', 'qualityscore', 'dfg', 'ac_helix',
                           'missing_residues', 'pocket']]

    # add column with positions of missing residues (replacing column with number of missing residues)
    df_screen = add_missing_residues(df_screen)

    # For each kinase with x different pdb codes: for each pdb code keep the structure with the best quality score
    df_screened = df_screen.groupby(["kinase", "pdb"]).max()["qualityscore"].reset_index()

    # Merge with df_screen, because we need chain information for next filtering step
    # left merge, as meta data from right frame is added to the left frame
    df_screened = pd.merge(df_screened, df_screen, how='left', on=['kinase', 'pdb', 'qualityscore'])

    # if same pdb code with same best score but different chain number: keep first chain
    df_screened.drop_duplicates(subset=["kinase", "pdb", "qualityscore"], keep='first', inplace=True)
    df_screened.reset_index(drop=True, inplace=True)

    # ------------ CHECKING --------------

    # # number of distinct pdb codes existing for each kinase (just for checking with database)
    # df_group = df_screened.groupby(["kinase"]).pdb.nunique().reset_index()
    # # print(df_group)
    # # check if there are still multiple pdb codes per kinase
    # if df_screened.shape[0] != df_group.pdb.sum():
    #     print('ERROR: Something went wrong. Multiple PDB codes per kinase!')
    #     sys.exit()
    # # Check if there is a unique pdb code for every kinase structure in the screened data set
    # df_group = df_screened.groupby(["kinase"]).pdb.nunique().reset_index()
    # if df_screened.shape[0] != df_group.pdb.sum():
    #     print('PDB codes not unique amongst kinases!')

    # # plot gap rate
    # missingResidues = []
    # for entry in df_screened.missing_residues:
    #     missingResidues.extend(entry)
    # plt.hist(missingResidues, density=True, bins=85, rwidth=0.8)
    # plt.title('Gap rate of 85 binding pocket residues')
    # plt.xticks(range(5, 85, 5))
    # plt.xlabel('KLIFS residue number')
    # plt.show()
    # sys.exit()

    # -------------------------------------

    return df_screened


def add_missing_residues(df):

    """

    Given a dataframe with a 'pocket' column listing the 85 residues of a kinase binding pocket,
    create a list of missing residues and add/replace the 'missing_residue' column

    Parameters
    ----------
    df: Pandas DataFrame
        input dataframe including a 'pocket' column

    Returns
    -------
    df: Pandas Dataframe
        output dataframe with added 'missing_residue' column

    """

    missing_residues = []
    for ix, row in df.iterrows():
        missing_residues.append(find_missing_residues(row.pocket))
    df['missing_residues'] = missing_residues
    return df


def find_missing_residues(sequence):  # , numMissing):

    """

    Given a sequence of binding pocket residues, find positions of missing residues ('_')
    (numMissing was used to save time, but this value is not always correct.)

    Parameters
    ----------
    sequence: Str
        85 binding pocket residues

    Returns
    -------
    missing_residues: list(int)
        list of missing residue indices

    """

    # iterate over sequence string
    missing_residues = []
    for i, aa in enumerate(sequence):
        # all missing residues found
        # if numMissing == 0:
        #     return missingResidues
        # missing residue
        if aa == '_':
            missing_residues.append(i+1)
            # numMissing -= 1
    return missing_residues


def fix_residue_numbers(pocket_mol2, missing_residues):

    """

    Fix residue IDs according to KLIFS in Mol2 instance of binding pocket
    (In the Mol2 file, the numbering was consecutive w/o considering missing residues)

    Parameters
    ----------
    pocket_mol2: Pandas DataFrame
        atom block from mol2 file representing the pocket (read using PandasMol2().read_mol2())
    missing_residues: list(int)
        positions (KLIFS numbering) of missing residues

    Returns
    -------
    pocket_mol2: Pandas DataFrame
        updated mol2 block

    """

    for res in missing_residues:
        pocket_mol2['res_id'].mask(pocket_mol2['res_id'] >= res, pocket_mol2['res_id'] + 1, inplace=True)
    return pocket_mol2


def get_folder_name(df):

    """

    Given information on a specific KLIFS structure, get relative path from KLIFS_download/ to corresponding files

    Parameters
    ----------
    df: Pandas Series
        dataframe with species, kinase, pdb, alt, and chain

    Returns
    -------
    folder: Str
        path from KLIFS_download/ to respective files

    """

    if df.alt == ' ':
        folder = df.species.upper()+'/'+df.kinase+'/'+df.pdb+'_chain'+df.chain
    else:
        folder = df.species.upper()+'/'+df.kinase+'/'+df.pdb+'_alt'+df.alt+'_chain'+df.chain

    return folder


def get_file_name(df):

    """

    Given information on a specific KLIFS structure, get output file name (for drawings)

    Parameters
    ----------
    df: Pandas Series
        dataframe with pdb, alt, and chain

    Returns
    -------
    file: Str
        output file name for this structure

    """

    if df.alt == ' ':
        file = df.pdb+'_chain'+df.chain
    else:
        file = df.pdb+'_alt'+df.alt+'_chain'+df.chain

    return file
