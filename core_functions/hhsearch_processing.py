
import pandas as pd
import multiprocessing


# take a list of strings and return counts of words separated by spaces
# ignores anything contained in regex blacklist expression
def calculate_label_counts(labels, blacklist='(protein)'):
    words = [label.split() for label in labels]
    words = [item for sublist in words for item in sublist]

    # remove common phrases from filter
    words = pd.Series(words)[~(pd.Series(words).str.contains(blacklist, regex=True))]

    word_counts = words.value_counts()

    return word_counts


# basic try to ammend the cropping of profile headers employed by hhsuite, capped at 138 or 142 chars??
# requires searchDF from parse_HHsuite or equivalent Query, Target dataframe
# requires global reference header mapping for both query and targets containing acc and header info as index

def create_hhsuite_header_mapping(searchDF, global_header_mapping):
    accs = []

    # append Target and Query columns
    entries = pd.concat([searchDF.Query, searchDF.Target]).unique()

    # iterate over entries and try to find a accession
    for hit in entries:
        # initial attempt by direct matching
        try:
            hit_acc = global_header_mapping.loc[hit].acc

        # for cropped entriestry to slice acc from first space separated element
        # then refer to the global mapping for header
        except KeyError:
            print('cannot find hit for', hit)
            print('trying via acc')
            hit_acc = hit.split()[0]
            new_header = global_header_mapping[global_header_mapping.acc == hit_acc].index[0]
            # print(new_header)
            # print('new acc is', hit_acc)
        accs.append(hit_acc)

    # format and return DF
    hhsuite_header_mapping = pd.DataFrame({'acc': accs, 'header': entries})

    hhsuite_header_mapping.sort_values(by='header', inplace=True)
    hhsuite_header_mapping.set_index(keys=['header'], drop=True, inplace=True)

    return hhsuite_header_mapping


# parse HHSuite outfut fromm ffdata into pkl files
def parse_and_write(file):
    print(file)
    new_data = HH.load_HHBlitsData(file)
    new_data.write_pkl(file + '.pkl')
    # new_data.write_data_tsv(file+'.tsv')

def parse_filter_write(file):
    thread = multiprocessing.current_process().pid
    print(f'{thread} reading {file}\n')
    new_data = HH.load_HHBlitsData(file)
    new_data.write_pkl(file + '.pkl')
    print(f'{thread} parsing\n')

    for key, query in new_data.data.items():
        query.add_self_hit()
        query.filter_numeric(field='Pairwise_cov', min=20, replace=True, keep_self=True)
        query.filter_numeric(field='Prob', min=50, replace=True, keep_self=True)

    print(f'{thread} writing\n')
    new_data.write_pkl(file + '.pkl_filtered')