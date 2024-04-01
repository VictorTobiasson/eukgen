import math
import pickle
import time

import numpy as np


#helper functions involved in IO
#small parsing tasks, general functions or repeated calculations

# SEQUENCE ALIGNMENTS OPERATIONS


# helper IO functions for basic fasta reading into {id:seq} dict
# id taken as string up until first space character
def fasta_to_dict(fastastring=None, file=None):
    if file != None:
        with open(file, 'r') as fastafile:
            fastalines = fastafile.readlines()

    # trim everything after first space in line. Avoids pathologic cases of headers such as "NR_XXXX (abc-->cd)"
    for n, line in enumerate(fastalines):
        fastalines[n] = ''.join(line.split(' ')[0])

    fastastring = '\n'.join(fastalines)

    entries = [entry.strip() for entry in ''.join(fastastring).split('>') if entry != '']
    fastas = {entry.split('\n')[0]: ''.join(entry.split('\n')[1:]) for entry in entries}

    # replace unknown characters with A
    # old blacklist: blacklist = set('BJUXZ*')
    blacklist = [';', '~', ':', 'u', '8', '[', '%', 'x', '^', '\x0c', '_',
                 ',', ' ', 'b', '\x0b', '.', '9', 'J', '1', '@', '}', '>',
                 '|', '\n', "'", '2', '$', 'z', '3', 'o', '<', '6', '=',
                 '#', '7', '{', 'Z', '&', '+', '(', '?', '/', '!', 'X',
                 'O', 'U', '0', '5', ')', ']', '*', 'B', '"', '4', '`', '\r', 'j', '\t']

    for key, seq in fastas.items():

        for char in set(seq):
            if char in blacklist:
                seq = seq.replace(char, 'A')
                #print(f'WARNING: {key} Replaced illegal {char} with A')

        fastas[key] = seq
    return fastas


# returns simple single line fasta from {id:seq} dict
def dict_to_fasta(seq_dict, write_file=False, verbose=True):
    fasta_str = '\n'.join(f'>{key}\n{value}' for key, value in seq_dict.items())

    if write_file != False:
        with open(write_file, 'w') as outfile:
            outfile.write(fasta_str)
        if verbose:
            print(f'Wrote {write_file}')

    return fasta_str


# small helper for pkl parsing
def load_pkl(pkl_file):
    with open(pkl_file, 'rb') as infile:
        item = pickle.load(infile)
    return item


def dump_pkl(item, pkl_file):
    with open(pkl_file, 'wb') as outfile:
        pickle.dump(item, outfile)
    print(f'Pickled item as {pkl_file}')

# pandas helper function to reset_index inplace
def reindex(df, column):
    df.sort_values(by=column, inplace=True)
    df.set_index(keys=[column], drop=True, inplace=True)


# define the entropy for a string given amino acids, protein=True or, DNA protein=False
def column_entropy(string, protein=True, gaptoken='-'):
    size = len(string)
    counts = [string.count(i) for i in set(string).difference({gaptoken})]
    entropy = -sum([i / size * math.log2(i / size) for i in counts])

    if protein:
        entropy_uniform = math.log2(20)
    else:
        entropy_uniform = 2

    gap_entropy = entropy_uniform * (string.count(gaptoken) / size)
    information = entropy_uniform - entropy - gap_entropy

    return max(information, 0)



# columnwise cut based on criteria
def filter_by_entropy(seq_dict, entropy_min, seq_length_frac_min=None, filter_accs=[], gaptoken='-'):
    # transpose seqs into columns
    cols = [''.join(seq) for seq in list(zip(*seq_dict.values()))]

    # filter only by columns present in filter_accs keys
    if filter_accs:
        filter_dict = {key: value for key, value in seq_dict.items() if key in filter_accs}
        filter_cols = [''.join(seq) for seq in list(zip(*filter_dict.values()))]

    else:
        filter_cols = cols

    # include based on entropy threshold
    filter_cols = [col for i, col in enumerate(cols) if column_entropy(filter_cols[i], gaptoken=gaptoken) > entropy_min]

    # transpose back to alignment
    filter_aln = [''.join(col) for col in list(zip(*filter_cols))]

    seq_dict = {key: value for key, value in zip(seq_dict.keys(), filter_aln)}

    # optional short sequence exclusion by fraction of total length
    if seq_length_frac_min is not None:
        max_len = len(list(seq_dict.values())[0])
        seq_dict = {key: value for key, value in seq_dict.items() if
                    len(value.replace(gaptoken, '')) >= max_len * seq_length_frac_min}

    # return with original keys
    return seq_dict


# input helper to read tsv from file, merge singletons
def read_cluster_tsv(cluster_file, split_large=False, max_size=500, batch_single=False, single_cutoff=1):
    # read TSV and group clusters based on first tsv column
    with open(cluster_file, 'r') as infile:
        clusters = {}

        for l in infile.readlines():
            cluster_acc, acc = l.strip().split('\t')

            if cluster_acc not in clusters.keys():
                clusters[cluster_acc] = [acc]

            else:
                clusters[cluster_acc].append(acc)

    # merge all clusters smaller than cutoff into one
    if batch_single:
        filter_dict = {}
        singles = []
        for key, accs in clusters.items():
            # gather singletons
            if len(accs) <= single_cutoff:
                singles.extend(accs)

            # keep larger clusters
            else:
                filter_dict[key] = accs

        if singles:
            filter_dict[singles[0]] = singles

        clusters = filter_dict

    # split clusters larger than x into smaller pieces
    if split_large:
        filter_dict = {}
        for key, accs in clusters.items():
            if len(accs) > max_size:
                # partition large clusters into batches of max_size
                for split in range(0, len(accs), max_size):
                    batch = accs[split:split + max_size]
                    filter_dict[batch[0]] = batch

            # keep smaller clusters
            else:
                filter_dict[key] = accs

        clusters = filter_dict

    return clusters


#calculate the entropy of a list of labels
def calculate_label_entropy(l):
    H = 0
    size = len(l)
    for i in set(l):
        n = l.count(i)
        hl = (n/size)*np.log2(1/(n/size))
        H += hl
    return H


#slurm run control
# get list of jobid statuses, every refresh cycle, return once all COMPLETED and none RUNNING or PENDING
def swarm_submit_with_lock(swarmfile, refresh=60, verbose=True):

    import subprocess
    from time import sleep, gmtime, strftime

    #run the swarmfile
    jobID = subprocess.run(f'swarm {swarmfile}'.split(), capture_output=True, text=True).stdout.strip()
    print(f'Swarm: {jobID} Submitted')

    #allow the sbatch system some time before checking
    running = True
    sleep(10)

    # wait for all swarms to finish
    while running:
        # timeout to not overburden the sbatch system
        sleep(refresh)
        # get list of all jobid statuses
        sacct_status = f'sacct -n -X -j {jobID} -o state%20 | sort | uniq -c'
        run_status = subprocess.run(sacct_status, capture_output=True, text=True, shell=True).stdout.strip()
        if verbose:
            t = time.gmtime()
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S", t)
            print(f'Swarm: {jobID} | {timestamp}\n{run_status}')

        if run_status.count('RUNNING') == 0 and run_status.count('COMPLETED') > 0 and run_status.count('PENDING') == 0:
            running = False

    if verbose:
        print(f'All statuses for {jobID} are "COMPLETED", returning.')

    return

#takes a list of tuples of (name: weight) and greedily pack them into bins
# bin capacity is either equal to the largest weight or total weight divided by number of bins with some margin
#returns list of item names and list of item weights
def greedy_pack(items, max_bins=1000):

    class Bin():
        #contains items as [(i, w_i), (j, w_j), ... (n, w_n)]
        def __init__(self):
            self.items = []
            self.weight = 0

        def add(self, item):
            self.items.append(item)
            self.weight = self.weight + item[1]

    #sort all items for faster packing
    items = sorted(items, key= lambda x: x[1])
    total_weight = sum(i[1] for i in items)

    #bin capacity is the largest element unless this causes to many bins
    capacity = items[-1][1]
    total_bins = int(total_weight / capacity)


    if total_bins > max_bins:
        print(f'With a bin capacity={capacity} given by the largest item there are more than {max_bins} bins ({total_bins})')
        capacity = int(total_weight / max_bins * 1.2)
        print(f'Estimating bin capacity total weight divided by {max_bins} (+20%) as {capacity}')

    #greedily pack bins
    packed_bins = [Bin()]

    while items:
        item = items.pop()
        #try to fit in all bins, break of found
        for bin in packed_bins:
            added = 0
            #add if space
            if bin.weight + item[1] <= capacity:
                bin.add(item)
                added = 1
                break

        #if never added make a new bin and add
        if added == 0:
            packed_bins.append(Bin())
            packed_bins[-1].add(item)

    packed_items = [[item[0] for item in bin.items] for bin in packed_bins]
    packed_weights = [[item[1] for item in bin.items] for bin in packed_bins]

    return packed_items, packed_weights










