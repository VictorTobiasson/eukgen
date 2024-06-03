from core_functions.mmseqs_functions import collapse_nonredundant, create_initial_profile_cluster, \
    cascade_profile_cluster
from core_functions.mmseqs_functions import extract_query_fasta, realign_all_fastas_swarm_binned
from core_functions.hhsuite_functions import build_HHsuite_database_from_MSAs, hhsuite_swarm_search

from paths_and_parameters import prok2111_as, euk72_ep, path_root, path_processing
from paths_and_parameters import realignment_swarm_opts, hhblits_swarm_opts, hhblits_search_opts

# # perform initial cascade clustering
root = '/data/tobiassonva/data/eukgen/processing/mini_euk_checked_sensitive/'
seqDB = '/data/tobiassonva/data/eukgen/core_data/mini_euk/mini_euk'
threads = 24

mmseqs_cascade_opts = {'cascade_steps': 4,
                       'nonredundant_cluster': '-s 2 -c 0.8 --cov-mode 0 --min-seq-id 0.9 --cluster-reassign 1 --max-iterations 1 --max-seqs 400 --remove-tmp-files',
                       'initial_cluster': '-s 7.5 -c 0.8 --cov-mode 0 -e 1e-10 --cluster-steps 2 --cluster-reassign 1 --max-seqs 1000 --remove-tmp-files',
                       'initial_search': '-s 7.5 -c 0.8 --cov-mode 0 --add-self-matches 1 -a --max-seqs 800 --num-iterations 3 --remove-tmp-files',
                       'result2profile': '-e 1E-10 --e-profile 1E-10 --cov 0.8',
                       'profile2consensus': '',
                       'cascade_search': '-s 7.5 -c 0.8 --cov-mode 0 --add-self-matches 1 -a --max-seqs 800 --num-iterations 2 --remove-tmp-files',
                       'cascade_clust': ''}

exe_realignment = path_root + 'muscle_crop_and_align.py'
muscle_realignment_timeout = 7200
realignment_swarm_opts = {'threads-per-process': 8,
                          'gb-per-process': '25',
                          'time': '72:00:00',
                          'logdir': '/data/tobiassonva/data/log/swarm_out/',
                          'job-name': 'swarm_realignment',
                          'maxrunning': 1000,

                          }

cascade_steps = mmseqs_cascade_opts['cascade_steps']


# make output file if not exists
import subprocess
subprocess.run(f'mkdir {root}'.split())

# reduce dataset be collapsing to 90% identity and 80% pairwise coverage
nonredundant_repeq = collapse_nonredundant(root, seqDB, threads, extra_opts=mmseqs_cascade_opts)

# initial cluster DB to form profiles
profile_clusterDB = create_initial_profile_cluster(root, nonredundant_repeq, threads, extra_opts=mmseqs_cascade_opts)

# cascade cluster DB
cascade_clusterDB = cascade_profile_cluster(root, nonredundant_repeq, profile_clusterDB, cascade_steps, threads,
                                            extra_opts=mmseqs_cascade_opts)

# extract fastas from cluters and realign
fasta_root = extract_query_fasta(root, seqDB, cascade_clusterDB, skip_singletons=True)

# realign all fastas using knapsacked swarm
fasta_root = realign_all_fastas_swarm_binned(fasta_root, realignment_swarm_opts, max_swarms=1000, max_sequences=100, filter_entropy=0.25, muscle_reps=1)

# # search new databse against itself
# output_DB_root = root
# HHsuiteDB_name = 'prok2311'
# mpis = 24
# search_root = root + 'self_search/'
# splits = 500
#
# fasta_root = '/vf/users/tobiassonva/data/prok2311/processing/cluster_fastas/'
#
# HHsuiteDB = build_HHsuite_database_from_MSAs(fasta_root, output_DB_root, HHsuiteDB_name, mpis=mpis)

# HHsuiteDB = '/vf/users/tobiassonva/data/eukgen/processing/euk72_ep3/euk72_ep'

# search_root = hhsuite_swarm_search(HHsuiteDB, HHsuiteDB, search_root, splits, hhblits_swarm_opts,
#                                    threads=mpis, search_type='hhblits', search_opts=hhblits_search_opts,
#                                    run_with_lscratch=True)
