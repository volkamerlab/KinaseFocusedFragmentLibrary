import pandas as pd


# code taken and adapted from
# https://github.com/AndreaVolkamer/KinaseSimilarity/blob/master/Bachelorarbeit/structFP_phchem/structFP_phchem_preprocessing.ipynb

def preprocess_klifs_data(path_to_klifs_download, path_to_klifs_export):

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


# replace number of missing residues with list of missing residues
def add_missing_residues(df):
    missing_residues = []
    for ix, row in df.iterrows():
        missing_residues.append(find_missing_residues(row.pocket))
    df['missing_residues'] = missing_residues
    return df


# find positions of missing residues
# numMissing was used to save time, but this value is not always correct. :(
def find_missing_residues(sequence):  # , numMissing):
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


# fix residue IDs according to KLIFS in mol2 instance of binding pocket
# missing_residues: positions (KLIFS numbering) of missing residues
def fix_residue_numbers(pocket_mol2, missing_residues):
    for res in missing_residues:
        pocket_mol2['res_id'].mask(pocket_mol2['res_id'] >= res, pocket_mol2['res_id'] + 1, inplace=True)
    return pocket_mol2


# input: 1D data frame with species, kinase, pdb, alt, and chain
# output: path from KLIFS download to respective files
def get_folder_name(df):

    if df.alt == ' ':
        folder = df.species.upper()+'/'+df.kinase+'/'+df.pdb+'_chain'+df.chain
    else:
        folder = df.species.upper()+'/'+df.kinase+'/'+df.pdb+'_alt'+df.alt+'_chain'+df.chain

    return folder


def get_file_name(df):

    if df.alt == ' ':
        file = df.pdb+'_chain'+df.chain
    else:
        file = df.pdb+'_alt'+df.alt+'_chain'+df.chain

    return file
