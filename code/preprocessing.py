import pandas as pd
import sys


# code taken and adapted from
# https://github.com/AndreaVolkamer/KinaseSimilarity/blob/master/Bachelorarbeit/structFP_phchem/structFP_phchem_preprocessing.ipynb

def preprocessKLIFSData(path_to_KLIFS_download, path_to_KLIFS_export):

    # read overview file
    df = pd.read_csv(path_to_KLIFS_download)
    # read export file
    df_csv = pd.read_csv(path_to_KLIFS_export)
    # rename columns
    df_csv.columns = ['kinase', 'family', 'groups', 'pdb', 'chain', 'alt', 'species', 'ligand', 'pdb_id',
                      'allosteric_name', 'allosteric_PDB', 'dfg', 'ac_helix']

    # # delete entries that are not present in both data frames
    # not_in_struc_download = df_csv[~df_csv['pdb'].isin(df['pdb'])]
    # while not_in_struc_download.shape[0] > 0:
    #     df_csv.drop(df_csv.index[[not_in_struc_download.index[0]]], inplace=True)
    #     df_csv.reset_index(drop=True, inplace=True)
    #     not_in_struc_download = df_csv[~df_csv['pdb'].isin(df['pdb'])]
    #     print('Deleted an entry in export table.')
    # not_in_csv_download = df[~df['pdb'].isin(df_csv['pdb'])]
    # while not_in_csv_download.shape[0] > 0:
    #     df.drop(df.index[[not_in_csv_download.index[0]]], inplace=True)
    #     df.reset_index(drop=True, inplace=True)
    #     not_in_csv_download = df[~df['pdb'].isin(df_csv['pdb'])]
    #     print('Deleted an entry in overview table.')

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
    df_screen = df_screen[['kinase', 'groups', 'species', 'pdb', 'pdb_id', 'alt', 'chain', 'qualityscore', 'dfg', 'ac_helix',
                           'missing_residues', 'pocket']]

    # add column with positions of missing residues (replacing column with number of missing residues
    missingResidues = []
    for ix, row in df_screen.iterrows():
        missingResidues.append(findMissingResidues(row.pocket, row.missing_residues))
    df_screen['missing_residues'] = missingResidues

    # For each kinase with x different pdb codes: for each pdb code keep the structure with the best quality score
    df_screened = df_screen.groupby(["kinase", "pdb"]).max()["qualityscore"].reset_index()

    # Merge with df_screen, because we need chain information for next filtering step
    # left merge, as meta data from right frame is added to the left frame
    df_screened = pd.merge(df_screened, df_screen, how='left', on=['kinase', 'pdb', 'qualityscore'])

    # if same pdb code with same best score but different chain number: keep first chain
    df_screened.drop_duplicates(subset=["kinase", "pdb", "qualityscore"], keep='first', inplace=True)
    df_screened.reset_index(drop=True, inplace=True)

    # ------------ CHECKING --------------

    # number of distinct pdb codes existing for each kinase (just for checking with database)
    df_group = df_screened.groupby(["kinase"]).pdb.nunique().reset_index()
    # print(df_group)
    # check if there are still multiple pdb codes per kinase
    if df_screened.shape[0] != df_group.pdb.sum():
        print('ERROR: Something went wrong. Multiple PDB codes per kinase!')
        sys.exit()
    # Check if there is a unique pdb code for every kinase structure in the screened data set
    df_group = df_screened.groupby(["kinase"]).pdb.nunique().reset_index()
    if df_screened.shape[0] != df_group.pdb.sum():
        print('PDB codes not unique amongst kinases!')

    # -------------------------------------

    return df_screened


# find positions of missing residues
def findMissingResidues(sequence, numMissing):
    # iterate over sequence string
    missingResidues = []
    for i, aa in enumerate(sequence):
        # all missing residues found
        if numMissing == 0:
            return missingResidues
        # missing residue
        if aa == '_':
            missingResidues.append(i+1)
            numMissing -= 1
    return missingResidues


def getFolderName(df):

    if df.alt == ' ':
        folder = df.species.upper()+'/'+df.kinase+'/'+df.pdb+'_chain'+df.chain
    else:
        folder = df.species.upper()+'/'+df.kinase+'/'+df.pdb+'_alt'+df.alt+'_chain'+df.chain

    return folder
