from core_functions.mmseqs_functions import collapse_nonredundant, create_initial_profile_cluster, cascade_profile_cluster
from core_functions.mmseqs_functions import extract_query_fasta, realign_all_fastas_swarm_binned
from core_functions.hhsuite_functions import build_HHsuite_database_from_MSAs, hhsuite_swarm_search

from paths_and_parameters import asgard2023, mmseqs_cascade_opts, path_processing
from paths_and_parameters import realignment_swarm_opts, hhblits_swarm_opts, hhblits_search_opts

# perform initial cascade clustering
root = path_processing + 'asgard2023/'
seqDB = asgard2023
threads = 48
cascade_steps = mmseqs_cascade_opts['cascade_steps']

# # make output file if not exists
# import subprocess
# subprocess.run(f'mkdir {root}'.split())

# # reduce dataset be collapsing to 90% identity and 80% pairwise coverage
# nonredundant_repeq = collapse_nonredundant(root, seqDB, threads, extra_opts=mmseqs_cascade_opts)

# # initial cluster DB to form profiles
# profile_clusterDB = create_initial_profile_cluster(root, nonredundant_repeq, threads, extra_opts=mmseqs_cascade_opts)

# # cascade cluster DB
# cascade_clusterDB = cascade_profile_cluster(root, nonredundant_repeq, profile_clusterDB, cascade_steps, threads,
#                                             extra_opts=mmseqs_cascade_opts)

cascade_clusterDB = '/data/luojaa/eukgen/processing/asgard2023/asgard2023_clean.repseq.cascaded_cluster'
# extract fastas from cluters and realign
# fasta_root = extract_query_fasta(root, seqDB, cascade_clusterDB)

# # reduce required threads as all subeseqent jobs are swarmed and specify own threads
# fasta_root = realign_all_fastas_swarm_binned(fasta_root, realignment_swarm_opts,
#                                              max_swarms=1000,
#                                              max_sequences=200,
#                                              filter_entropy=0.25,
#                                              muscle_reps=5)

# search new databse against itself
fasta_root = root + 'cluster_fastas/'
output_DB_root = root
HHsuiteDB_name = 'asgard2023_clean'
mpis = 24
search_root = '/data/luojaa/kegg_hmm_db/search_asgard/'
splits = 500

# HHsuiteDB = build_HHsuite_database_from_MSAs(fasta_root, output_DB_root, HHsuiteDB_name, mpis=mpis)

HHsuiteQueryDB = '/data/luojaa/eukgen/processing/asgard2023/asgard2023_clean'
HHsuiteTargetDB = '/data/luojaa/eukgen/processing/prok2111/prok2111'

search_root = hhsuite_swarm_search(HHsuiteQueryDB, HHsuiteTargetDB, search_root, splits, hhblits_swarm_opts,
                                   threads=mpis, search_type='hhblits', search_opts=hhblits_search_opts,
                                   run_with_lscratch=True)
