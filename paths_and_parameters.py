# main parameter storage

# global paths
path_root = '/data/tobiassonva/data/eukgen/'
path_core_data = path_root + 'core_data/'

path_processing = path_root + 'processing/'
path_tmp = path_root + 'tmp/'
path_testing = path_root + 'testing/'
path_notebooks = path_root + 'notebooks/'

# mmseqs_scratch
path_mmseqs_tmp = path_tmp + 'mmseqs_scratch/'

# databases

# core mmseqs euk72
path_euk72 = path_core_data + 'euk72/'
euk72 = path_euk72 + 'euk72'
euk72_filtered = path_euk72 + 'euk72_filtered'
euk72_header = path_euk72 + 'euk72_header_mapping.pkl'

# core mmseqs prok2111
path_prok2111 = path_core_data + 'prok2111/'
prok2111 = path_prok2111 + 'prok2111'
prok2111_filtered = path_prok2111 + 'prok2111_filtered'
prok2111_header = path_prok2111 + 'prok2111_header_mapping.pkl'

# prok2111 with asgards from ettema2023
path_prok2111_as = path_core_data + 'prok2111_as/'
prok2111_as = path_prok2111_as + 'prok2111_as'
prok2111_as_tax = path_prok2111_as + 'prok2111_as.tax'

# euk72 merged with cleaned and parsed eukprot
# already filtered 20 < len(seq) < 5000
path_euk72_ep = path_core_data + 'euk72_ep/'
euk72_ep = path_euk72_ep + 'euk72_ep'
euk72_ep_tax = path_euk72_ep + 'euk72_ep.tax'



# small prok dataset for testing ~9000 proteins
path_prok_small = path_core_data + 'small_prok2111_subset/small_prok'

# taxonomy
path_taxonomy = path_core_data + 'taxonomy/'
ncbi_taxonomy = path_taxonomy + 'ncbi_tax_231127/'

# merged tax files from prok2111_as and euk72_ep with columns ['acc', 'orgid', 'superkingdom', 'class']
merged_protein_tree_taxonomy = path_taxonomy + 'euk_prok_merged_protein.tax'

# protein tax mapping and species lienage mapping from NCBI
euk72_protein_taxonomy = path_taxonomy + 'euk72_protein_taxonomy.pkl'
euk72_species_taxonomy = path_taxonomy + 'euk72_species_taxa.tsv'
prok2111_protein_taxonomy = path_taxonomy + 'prok2111_protein_taxonomy.pkl'
prok2111_species_taxonomy = path_taxonomy + 'prok2111_species_taxa.tsv'

# smaller streamlined version only covering proteins found in euk searches from July
prok2111_protein_taxonomy_trimmed = path_taxonomy + 'prok2111_protein_taxonomy_trimmed.pkl'

# executable paths
exe_python = 'python'
exe_famsa = '/data/tobiassonva/data/software/FAMSA-2.0.1/famsa'
exe_fasttree = 'FastTree'
exe_muscle5 = 'muscle'
exe_iqtree = 'iqtree2'
exe_mmseqs = 'mmseqs'
exe_hhsearch_omp = 'hhsearch_omp'
exe_hhblits_omp = 'hhblits_omp'

# mmseqs cascaded clustering parameters
mmseqs_cascade_opts = {'cascade_steps': 4,
                       'nonredundant_cluster': '-s 2 -c 0.8 --cov-mode 0 --min-seq-id 0.9 --cluster-reassign 1 --max-iterations 1 --max-seqs 1000 --remove-tmp-files',
                       'initial_cluster': '-s 5 -c 0.8 --cov-mode 0 -e 1e-10 --cluster-steps 2 --cluster-reassign 1 --max-seqs 1000 --remove-tmp-files',
                       'initial_search': '-s 6 -c 0.8 --cov-mode 0 --add-self-matches 1 -a --max-seqs 300 --num-iterations 2 --remove-tmp-files',
                       'result2profile': '-e 1E-10 --e-profile 1E-10 --cov 0.8',
                       'profile2consensus': '',
                       'cascade_search': '-s 7.5 -c 0.8 --cov-mode 0 --add-self-matches 1 -a --max-seqs 300 --num-iterations 2 --remove-tmp-files',
                       'cascade_clust': ''
                       }

exe_realignment = path_root + 'muscle_crop_and_align.py'
muscle_realignment_timeout = 7200
realignment_swarm_opts = {'threads-per-process': 8,
                          'gb-per-process': '50',
                          'time': '72:00:00',
                          'logdir': '/data/tobiassonva/data/log/swarm_out/',
                          'job-name': 'swarm_realignment',
                          'maxrunning': 1000,

                          }

# hhsuite database formatting
hhconsensus_opts = '-M 50 -maxres 65535 -v 0'
# hhmake requires at least 2 shown sequences or else it segfaults...?
hhmake_opts = '-seq 2 -v 0'
cstranslate_opts = '-f -x 0.3 -c 4 -I a3m'

# hhsuite searches
exe_hhsuite_search = path_root + 'hhsuite_search.py'
#note that the alignment need to be present for parsing of full names using parse_HHsuite
hhblits_search_opts = '-n 1 -p 80 -z 1 -Z 500 -B 500 -b 1'
hhblits_swarm_opts = {'threads-per-process': 24,
                      'gb-per-process': '100',
                      'time': '48:00:00',
                      'gres': 'lscratch:100',
                      'logdir': '/home/tobiassonva/data/log/swarm_out/',
                      'job-name': 'hhsuite_search',
                      'maxrunning': 500,
                      }

#more lscratch memory to accomodate uniref
hhblits_swarm_opts_uniref = {'threads-per-process': 24,
                      'gb-per-process': '100',
                      'time': '48:00:00',
                      'gres': 'lscratch:300',
                      'logdir': '/data/tobiassonva/data/log/swarm_out/',
                      'job-name': 'hhsuite_search',
                      'maxrunning': 500,
                      }

#iqtree saturates at 8 threads for most alignments
microcosm_format_opts = {'original_query_DB': euk72_ep,
                         'original_target_DB': prok2111_as,
                         'taxonomy_mapping': merged_protein_tree_taxonomy,
                         'max_euk_sequences': 30,
                         'max_prok_sequences': 100,
                         'filter_entropy': 0.25,
                         'filter_length_frac': 0.2,
                         'threads': 8,
                         'muscle_reps': 1,
                         'muscle_timeout': 1800,
                         'evo_model_params': '-m MFP -mset LG,Q.pfam --cmin 4 --cmax 12'
                         }