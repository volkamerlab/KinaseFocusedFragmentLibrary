import pandas as pd


# code taken and adapted from
# https://github.com/AndreaVolkamer/KinaseSimilarity/blob/master/Bachelorarbeit/structFP_phchem/structFP_phchem_preprocessing.ipynb


def read_klifs_meta_data(path_to_klifs_download, path_to_klifs_export):

    """
    Given the two CSV files downloaded from KLIFS:
    - merges the two dataframes
    - selects relevant columns
    - adds a column listing the missing residue numbers

    Parameters
    ----------
    path_to_klifs_download: Path or Str
        path to overview.csv file given by KLIFS (including the file name)
    path_to_klifs_export: Path or Str
        path KLIFS_export.csv file (including the file name)

    Returns
    -------
    df: Pandas DataFrame
        Dataframe including merged KLIFS metadata

    """

    # read overview file
    df_overview = pd.read_csv(path_to_klifs_download)
    # read export file
    df_export = pd.read_csv(path_to_klifs_export)
    # rename columns
    df_export.columns = ['kinase', 'family', 'group', 'pdb', 'chain', 'alt', 'species', 'ligand', 'pdb_id',
                      'allosteric_name', 'allosteric_PDB', 'dfg', 'ac_helix']

    # sync kinase names for both data frames
    for ix, row in df_export.iterrows():
        alt = row.alt
        kinase = row.kinase
        if '(' in kinase:
            row.kinase = kinase[kinase.find("(") + 1:kinase.find(")")]
        if alt == '-':
            row.alt = ' '

    # outer merge, as information from both data frames should be kept
    df = pd.merge(df_overview, df_export, how='outer', on=['species', 'kinase', 'pdb', 'chain', 'alt', 'allosteric_PDB'])

    # loose irrelevant data
    df = df[['kinase', 'family', 'group', 'species', 'pdb', 'pdb_id', 'alt', 'chain', 'qualityscore', 'dfg', 'ac_helix',
                           'missing_residues', 'pocket']]

    # add column with positions of missing residues (replacing column with number of missing residues)
    df = add_missing_residues(df)

    return df


def choose_best_klifs_structure(klifs_meta_data):

    """
    Chooses the best quality KLIFS structure for each PDB in a dataframe containing metadata from KLIFS

    Parameters
    ----------
    klifs_meta_data: Pandas DataFrame
        DataFrame containing metadata from KLIFS, created using read_klifs_meta_data

    Returns
    -------
    klifs_meta_data_filtered: Pandas DataFrame
        Filtered dataframe including only the best quality structure for each PDB

    """

    # For each kinase with x different pdb codes: for each pdb code keep the structure with the best quality score
    klifs_meta_data_filtered = klifs_meta_data.groupby(["kinase", "pdb"]).max()["qualityscore"].reset_index()

    # Merge with df_screen, because we need chain information for next filtering step
    # left merge, as meta data from right frame is added to the left frame
    klifs_meta_data_filtered = pd.merge(klifs_meta_data_filtered, klifs_meta_data, how='left', on=['kinase', 'pdb', 'qualityscore'])

    # if same pdb code with same best score but different chain number: keep first chain
    klifs_meta_data_filtered.drop_duplicates(subset=["kinase", "pdb", "qualityscore"], keep='first', inplace=True)
    klifs_meta_data_filtered.reset_index(drop=True, inplace=True)

    return klifs_meta_data_filtered


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
