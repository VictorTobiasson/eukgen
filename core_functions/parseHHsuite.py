#!/usr/bin/env python
# coding: utf-8

# standalone python only module to reformat hhblits and hhsearch ffdata into a python data structure
# requires pickle and/or pandas for output

# main class for aggregates queries
class HHblitsData:
    def __init__(self):
        self.data = {}
        self.size = len(self.data.keys())
        self.query_names = list(self.data.keys())

    # add library of {name:HHquery, ...}
    def add_entries(self, HHblitsData_data):
        self.data.update(HHblitsData_data)
        self.size = len(self.data.keys())
        self.query_names = list(self.data.keys())

    def load_from_HHR(self, file):
        new_data = load_HHBlitsData(file)
        self.add_entries(new_data.data)

    # tricky cases as depends on what is in the pickle file?
    # currently load pickle files of whole HHblitsData objects
    def load_from_pkl(self, file):
        import pickle
        with open(file, 'rb') as infile:
            self.add_entries(pickle.load(infile).data)

    def write_data_tsv(self, filename):
        import os
        # write the header info
        header = 'Query\t' + '\t'.join([key for key in self.data[self.query_names[0]].hit_dict.keys()]) + '\n'

        # write each entry to file
        with open(filename, 'w') as outfile:
            outfile.write(header)

            for key, item in self.data.items():
                outfile.write(item.format_data_tsv() + '\n')
        cut_filename = filename.replace(".tsv", ".cut.tsv")
        os.system(f"cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 {filename} > {cut_filename}")
        os.system(f"rm {filename}")
                

    def write_query_tsv(self, filename):
        # write header info
        header = '\t'.join([str(key) for key in self.data[self.query_names[0]].query_dict.keys()])

        # write each entry to file
        with open(filename, 'w') as outfile:
            outfile.write(header + '\n')

            for key, item in self.data.items():
                outfile.write('\t'.join([str(value) for entry, value in item.query_dict.items()]) + '\n')

    def write_pkl(self, filename):
        import pickle
        with open(filename, 'wb') as outfile:
            pickle.dump(self, outfile)


# class for single qury run of hhblits (should parse with HHsearch as well)
class HHblitsQuery:
    def __init__(self, query_dict, hit_dict):
        self.query_dict = query_dict
        self.hit_dict = hit_dict
        self.name = self.query_dict['Query']
        self.size = len(hit_dict['Target'])

    def format_data_tsv(self):
        # for multiple hits
        if isinstance(self.hit_dict['Target'], list):
            # transform into wide form for printing
            hit_items_T = transpose_list_of_lists([value for key, value in self.hit_dict.items()])

            # tsv join and add query name as first column
            tsv_lines = '\n'.join([self.name + '\t' + '\t'.join([str(n) for n in item]) for item in hit_items_T])

        # for single hits
        else:
            tsv_lines = self.name + '\t' + '\t'.join([str(value) for key, value in self.hit_dict.items()])

        return tsv_lines

    # simple numeric filter for entry values
    def filter_numeric(self, field, min=None, max=None, replace=False, keep_self=False):

        # create a boolean list of filter conditions
        if min and max:
            filter_list = [True if min < value < max else False for value in self.hit_dict[field]]
        elif min:
            filter_list = [True if min < value else False for value in self.hit_dict[field]]
        else:
            filter_list = [True if value < max else False for value in self.hit_dict[field]]

        # save self hits from filter by setting true
        if keep_self:
            if self.name not in self.hit_dict['Target']:
                print(f'{self.name} has no self hit! Will not filter as keep_self=True')
                return None
            else:
                self_hit_index = [i for i, value in enumerate(self.hit_dict['Target']) if value == self.name]
                filter_list = [entry if i not in self_hit_index else True for i, entry in enumerate(filter_list)]

        # keep all entries at indices from rows which meet the condition
        new_hit_dict = {key: [n for i, n in enumerate(value) if filter_list[i]] for key, value in self.hit_dict.items()}

        # update in place or return
        if replace:
            self.hit_dict = new_hit_dict
            self.size = len(self.hit_dict['Target'])


        else:
            return new_hit_dict

    def filter_pairwise_coverage(self, cutoff, replace=False, keep_self=False):
        # pairwise minimum filter
        pairwise_min_cov = calculate_pairwise_coverage(self.hit_dict)
        filter_list = [True if value > cutoff else False for value in pairwise_min_cov]

        # save self hits from filter by setting true
        if keep_self:
            if self.name not in self.hit_dict['Target']:
                print(f'{self.name} has no self hit! Will not filter as keep_self=True')
                return None
            else:
                self_hit_index = [i for i, value in enumerate(self.hit_dict['Target']) if value == self.name]
                filter_list = [entry if i not in self_hit_index else True for i, entry in enumerate(filter_list)]

        # keep all entries at indices from rows which meet the condition
        new_hit_dict = {key: [n for i, n in enumerate(value) if filter_list[i]] for key, value in self.hit_dict.items()}

        # update in place or return
        if replace:
            self.hit_dict = new_hit_dict
            self.size = len(self.hit_dict['Target'])
        else:
            return new_hit_dict

    # add mock self hits to the query with fake estimates of statistics
    def add_self_hit(self, force=False):
        # if query name among target names
        if self.name in self.hit_dict['Target'] and not force:
            return
            # print('Self hit exists')

        else:
            # print('Generating self hit')
            self_hit_values = ['FAKESELFHIT'] * len(self.hit_dict)
            query_len = self.query_dict['Match_columns']
            query_neff = self.query_dict['Neff']

            self_hit_values[0] = self.name
            self_hit_values[1:17] = [100.0, 0, 0, 99999, 0, query_len, 100, 100, 99999, 1, query_len, 1, query_len,
                                     query_len, query_neff, 1]

            self_hit = {entry[0]: self_hit_values[i] for i, entry in enumerate(self.hit_dict.items())}

            # merge the dicts again and update
            self.hit_dict = {key: [self_hit[key]] + value for key, value in self.hit_dict.items()}

    # TODO: finish styling
    def print_alignment(self, index):
        td = [value[index] for key, value in self.hit_dict.items()]
        print(self.name + '\n' + "\n".join(td[-6:]) + '\n' + td[0])


#helper function for creating a list of profile object from a single ffdata file
def loadHHsuiteProfiles(file):
    preformat_data = preformat_HHsuiteFFdata(file)
    profile_dicts = [parse_HHsuiteProfile(profile) for profile in preformat_data]
    profiles = [HHsuiteProfile(*dicts) for dicts in profile_dicts]
    return profiles


# main profile class
class HHsuiteProfile:
    def __init__(self, header_dict, consensus_dict, HMM_dict):
        self.hmm_dict = HMM_dict
        self.header_dict = header_dict
        self.consensus_sequence = consensus_dict['consensus_sequence']
        self.consensus_model = consensus_dict['consensus_model']
        self.name = self.header_dict['NAME']
        self.length = len(HMM_dict['A'])

    def write_data_tsv(self, filename):
        hmm_data = transpose_list_of_lists([value for value in self.hmm_dict.values()])
        header = '\t'.join(self.hmm_dict.keys())
        hmm_data = ['\t'.join(line) for line in hmm_data]
        tsv_lines = header + '\n' + '\n'.join(hmm_data)

        with open(filename, 'w') as outfile:
            outfile.write(tsv_lines)

    def write_pkl(self, filename):
        import pickle
        with open(filename, 'wb') as outfile:
            pickle.dump(self, outfile)

# helper function create an HHblitsData object from HHR file.
def load_HHBlitsData(file):
    # format list of queries
    preformat_data = preformat_HHsuiteFFdata(file)
    entry_dicts = [parse_HHsuiteHHR(entry) for entry in preformat_data]
    queries = {query_dict['Query']: HHblitsQuery(query_dict, hit_dict) for query_dict, hit_dict in entry_dicts}

    # add queries and return
    new_data = HHblitsData()
    new_data.add_entries(queries)

    return new_data

# helper function to convert [[a,b],[a,b],[a,b]] to [[a,a,a],[b,b,b]]
def transpose_list_of_lists(list_of_lists):
    return list(map(list, zip(*list_of_lists)))

# takes a HHsuite FFdata and preformats it for parsing with parse_HHsuiteHHR
def preformat_HHsuiteFFdata(file):
    with open(file) as ffdata:
        # read data as lines and merge
        queries = ''.join(ffdata.readlines())

        # search and replace for annoying words
        replace_dict = {'Query HMM': 'Query-HMM',
                        'Template HMM': 'Template-HMM'}

        for key, value in replace_dict.items():
            queries = queries.replace(key, value)

        # split at null delimiter to extract entries, drop last empty entry
        queries = queries.split('\x00')[0:-1]
        entries = [q.split('\n') for q in queries]

        # separate each entry into lines and drop empty lines
        entries = [[line for line in entry if line != ''] for entry in entries]

        return entries

# takes a HHsuite FFdata and preformats it for parsing with parse_HHsuiteHHR
def preformat_HHsuiteFFdata_Serial(file, entry_limit=99999999):

    # initalize values
    entries_read = 0
    entries = []
    entry = ''

    # search and replace for annoying words
    replace_dict = {'Query HMM': 'Query-HMM',
                    'Template HMM': 'Template-HMM'}

    with open(file, 'r') as ffdata:
        print('Opened', file)

        # read until limit
        while entries_read < entry_limit:
            line = ffdata.readline()

            # if line begins will null
            if line[0] == '\x00':

                print(f'Found entry number {entries_read+1}')

                # replace all bad names
                for key, value in replace_dict.items():
                    entry = entry.replace(key, value)

                # add entry to entries and start new entry
                entries.append([line for line in entry.split('\n') if line != ''])
                entries_read += 1
                entry = line.strip('\x00')

            else:
                entry += line

    return entries


# quite long procedural, could be heavily functionally refactored
def parse_HHsuiteHHR(entry):
    # PARSE THE QUERY BLOCK
    query_block = entry[0:7]
    query_data = [[n.split(' ')[0], ' '.join(n.split(' ')[1:]).strip()] for n in query_block]
    query_dict = {entry[0]: entry[1] for entry in query_data}

    # clean query_dict data formats
    for key, value in query_dict.items():
        if key in ['Match_columns', 'Neff', 'Searched_HMMs']:
            query_dict[key] = float(value)

    # PARSE THE HIT BLOCK

    # parse 7th line of labels and drop first 2 as entry number and cropped target label is skipped below
    hit_labels = [label.strip() for label in entry[7].split()][2:]
    hit_labels.append('Template_columns')

    hit_block = []
    hit_num = 0
    hit_block_T = [[] for _ in hit_labels]

    # iterate until first 3 characters is not a number
    for line in entry[8:]:
        try:
            hit_num = int(line[0:3])
            isint = True
        except ValueError:
            isint = False

        if isint:
            # parse hit line by skipping initial text block [0:35],
            # splitting, removing () and null entries and
            # turning to values to float
            line_trimmed = line[35:]
            line_split, template_columns = line_trimmed.split('(')
            data_list = [n for n in line_split.split(' ') if n != ''] + [template_columns[:-1]]

            hit_block.append(data_list)

            # transpose list of lists for dict
            hit_block_T = transpose_list_of_lists(hit_block)

        else:
            break

    # add all data to returned dict
    hit_dict = {label: hit_block_T[index] for index, label in enumerate(hit_labels)}

    # ITERATIVELY PARSE ALIGNMENT BLOCKS

    # get start of alignment blocks as where the line start with 'No '
    # add -1 as the final end value for later iteration
    alignment_starts = [line_no for line_no, line in enumerate(entry) if line[:3] == 'No '] + [-1]

    # get label names from second row after first entry
    # drop the first 4 as they are redundant with the data found in the hit table above
    target_labels = [value.split('=')[0] for value in entry[alignment_starts[0] + 2].split() if value != ''][4:]
    target_names = []
    target_data_list = []

    # get the initial alignment block width as length of third rightmost element
    # get the initial alignment block starting index as the match block length minus alignment block width
    # neccesary for cropping alignments to length dynamically
    first_query_line = entry[alignment_starts[0] + 3]
    first_match_line = entry[alignment_starts[0] + 5]
    alignment_width = len([value for value in first_query_line.split(' ') if value != ''][-3])
    alignment_start = len(first_match_line) - alignment_width

    # manually assign alignment data labels and initialize the array of alignment list types
    # for secondary structure assignments this would have to be added here together with some error handling
    alignment_labels = ['Query_sequence', 'Query_consensus', 'Matches', 'Target_consensus', 'Target_sequence',
                        'Confidence']
    alignment_data_list = [[], [], [], [], [], []]

    # for each alignment block start value, get the range of lines until the start value.
    # for each such start process the target_data and alignment_blocks
    for i in range(len(alignment_starts) - 1):
        alignment = entry[alignment_starts[i]:alignment_starts[i + 1]]

        # get the name of the target as the first line after "No X"
        target_names.append(alignment[1].strip('>'))

        # add the data from the second line skipping the first 4 redundant entries
        target_data_list.append([value.split('=')[-1] for value in alignment[2].split() if value != ''][4:])

        # iterate though alignment lines and concat all alignment
        alignment_data = ['', '', '', '', '', '']

        for line in alignment:
            # for flow control
            add_data = True

            # check for each type of recognized data line start
            if line[:3] == 'Q C':
                update_column = 1

            elif line[0] == 'Q':
                update_column = 0

                # lengths of alignments vary in the last section of each alignment block
                # check the length before the first line in each block
                alignment_width = len([value for value in line.split(' ') if value != ''][-3])

            elif line[0] == ' ':
                update_column = 2

            elif line[:3] == 'T C':
                update_column = 3

            elif line[0] == 'T':
                update_column = 4

            elif line[0] == 'C':
                update_column = 5

            else:
                # if a row is not recognized dont add any data
                add_data = False

            # cut the alignment information and extend the correct column
            if add_data:
                alignment_data[update_column] += line[alignment_start:alignment_start + alignment_width]

        # for every block append the corresponding columns to the alignment data
        alignment_data_list = [data + [alignment_data[i]] for i, data in enumerate(alignment_data_list)]

    # transpose the final target_data_list and add it to the hit_dict together with the names
    target_data_list_T = transpose_list_of_lists(target_data_list)
    hit_dict.update({label: target_data_list_T[index] for index, label in enumerate(target_labels)})
    hit_dict.update({'Target': target_names})

    # add the alignment data to the hit_dict
    hit_dict.update({label: alignment_data_list[index] for index, label in enumerate(alignment_labels)})

    # CLEANUP

    # reformat some messy data columns for easier use with Pandas
    extra_dict = {}
    final_key_order = ['Target', 'Prob', 'E-value', 'P-value', 'Score',
                       'SS', 'Cols', 'Identities', 'Similarity', 'Sum_probs',
                       'Query-HMM-start', 'Query-HMM-end', 'Template-HMM-start',
                       'Template-HMM-end', 'Template_columns', 'Template_Neff', 'Pairwise_cov', 'Query_sequence',
                       'Query_consensus', 'Matches', 'Target_consensus', 'Target_sequence', 'Confidence']

    # for entries without hits parsed (typically empty) return an empty library
    if hit_num == 0:
        return query_dict, {key: [] for key in final_key_order}

    for key, value in hit_dict.items():

        # convert scalar values to float or int
        if key in ['Prob', 'E-value', 'P-value', 'Score', 'SS', 'Similarity', 'Sum_probs', 'Template_Neff',
                   'Pairwise_cov']:
            hit_dict[key] = [float(n) for n in value]
        if key in ['Cols', 'Template_columns']:
            hit_dict[key] = [int(n) for n in value]

            # convert percentages
        if key == 'Identities':
            hit_dict[key] = [int(n.strip('%')) for n in value]

        # split column ranges and crop out HMM column from template HMM hybrid entry
        elif key in ['Query-HMM', 'Template-HMM']:
            data_list = [[int(a) for a in data.split('-')] for data in value]
            start, end = transpose_list_of_lists(data_list)
            extra_dict[key + '-start'] = start
            extra_dict[key + '-end'] = end

    # delete redundant columns and add new ones and reorder final dict
    hit_dict.pop('Query-HMM', None)
    hit_dict.pop('Template-HMM', None)
    hit_dict.update(extra_dict)

    # add pairwise coverage in place and reorder columns for display
    hit_dict['Pairwise_cov'] = calculate_pairwise_coverage(hit_dict, query_dict)

    hit_dict = {key: hit_dict[key] for key in final_key_order}

    # differnt amounts of hits and alignments can be printed by hhsuite which is a problem
    # pad all hit_dict data series to same lengths with 'nan' so it can be converted to DataFrame
    # remove one from alignment block number as it is padded with :-1 for iterating
    max_data_len = max(hit_num, len(alignment_starts) - 1)
    for key, value in hit_dict.items():
        hit_dict[key] = value + [float('nan')] * (max_data_len - len(value))

    return query_dict, hit_dict


def calculate_pairwise_coverage(hit_dict, query_dict):
    qs = hit_dict['Query-HMM-start']
    qe = hit_dict['Query-HMM-end']
    qn = query_dict['Match_columns']

    ts = hit_dict['Template-HMM-start']
    te = hit_dict['Template-HMM-end']
    tn = hit_dict['Template_columns']

    # (q_end-q_start)/t_columns
    ql = [i - j + 1 for i, j in zip(qe, qs)]
    qt_cov = [i / j for i, j in zip(ql, tn)]

    # can be calculated directly as query columns are same
    tq_cov = [(i - j + 1) / qn for i, j in zip(te, ts)]

    # pairwise minimum filter
    pairwise_min_cov = [min(i, j) for i, j in zip(qt_cov, tq_cov)]

    return pairwise_min_cov


# parse _hmm.ffdata part of DB for visualisation and statistics.
def parse_HHsuiteProfile(entry):
    header_stop = entry.index('SEQ')
    consensus_stop = entry.index('#')

    # parse header
    header_dict = {line.split()[0]: ' '.join(line.split()[1:]).strip() for line in entry[1:header_stop]}

    # parse consensus
    merged_consensus = [seq for seq in '\n'.join(entry[header_stop + 1:consensus_stop]).split('>') if seq != '']

    consensus_headers = ['consensus_model', 'consensus_sequence']
    consensus_dict = {consensus_headers[i]: ''.join(seq.split('\n')[1:]) for i, seq in enumerate(merged_consensus)}

    # parse state values
    NULL = [value.strip() for value in (entry[consensus_stop + 1].strip('NULL') + entry[consensus_stop + 4]).split('\t')
            if value != '']
    HMM = [value.strip() for value in (entry[consensus_stop + 2].strip('HMM') + entry[consensus_stop + 3]).split('\t')
           if value != '']
    HMM = [''.join(entry.split('->')) for entry in HMM]

    # add all data rows
    row_values = []

    for line in entry[consensus_stop + 5:-1]:
        if line[0] != ' ':
            row_values.append(line[7:].split('\t')[:-1])
        else:
            row_values[-1].extend(line[7:].split('\t')[:-1])

    # format dictionary with HMM header and rows
    HMM_dict = {HMM[i]: value for i, value in enumerate(transpose_list_of_lists(row_values))}

    return header_dict, consensus_dict, HMM_dict
