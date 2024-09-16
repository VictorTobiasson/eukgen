# main parameter storage

# global paths
path_root = '/vf/users/luojaa/eukgen/'
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

# core mmseqs asgard2023
path_asgard2023 = path_core_data + 'asgard2023_clean/'
asgard2023 = path_asgard2023 + 'asgard2023_clean'

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
merged_protein_tree_taxonomy =  '/data/luojaa/eukgen/mmseqs_victor/euk_prok_merged_protein_revised.tax'
#merged_protein_tree_taxonomy = "/data/luojaa/taxids/kegg_new_classes.mcrcsm.tsv"

# protein tax mapping and species lienage mapping from NCBI
euk72_protein_taxonomy = path_taxonomy + 'euk72_protein_taxonomy.pkl'
euk72_species_taxonomy = path_taxonomy + 'euk72_species_taxa.tsv'
prok2111_protein_taxonomy = path_taxonomy + 'prok2111_protein_taxonomy.pkl'
prok2111_species_taxonomy = path_taxonomy + 'prok2111_species_taxa.tsv'

# smaller streamlined version only covering proteins found in euk searches from July
prok2111_protein_taxonomy_trimmed = path_taxonomy + 'prok2111_protein_taxonomy_trimmed.pkl'

# executable paths
exe_python = 'python'
exe_famsa = '/data/luojaa/software/FAMSA/famsa'
exe_fasttree = 'FastTree'
exe_muscle5 = 'muscle'
exe_iqtree = 'iqtree2'
exe_mmseqs = 'mmseqs'
exe_hhsearch_omp = 'hhsearch_omp'
exe_hhblits_omp = 'hhblits_omp'
exe_foldseek = 'foldseek'
exe_icarus = 'python -u /home/tobiassonva/data/software/ICARUS/icarus.py'
exe_icarus_environment = 'icarus'

# mmseqs cascaded clustering parameters
mmseqs_cascade_opts = {'cascade_steps': 4,
                       'nonredundant_cluster': '-s 3 -c 0.8 --cov-mode 0 --min-seq-id 0.9 --cluster-reassign 1 --max-iterations 1 --max-seqs 1000 --remove-tmp-files',
                       'initial_cluster': '-s 4 -c 0.8 --cov-mode 0 -e 1e-10 --cluster-steps 2 --cluster-reassign 1 --max-seqs 1000 --remove-tmp-files',
                       'initial_search': '-s 5 -c 0.8 --cov-mode 0 --add-self-matches 1 -a --max-seqs 1000 --num-iterations 2 --remove-tmp-files',
                       'result2profile': '-e 1E-10 --e-profile 1E-10 --cov 0.8',
                       'profile2consensus': '',
                       'cascade_search': '-s 7.5 -c 0.8 --cov-mode 0 --add-self-matches 1 -a --max-seqs 500 --num-iterations 1 --remove-tmp-files',
                       'cascade_clust': ''
                       }

exe_realignment = path_root + 'muscle_crop_and_align.py'
muscle_realignment_timeout = 4000
realignment_swarm_opts = {'threads-per-process': 8,
                          'gb-per-process': '50',
                          'time': '72:00:00',
                          'logdir': '/data/luojaa/log/swarm_out/',
                          'job-name': 'swarm_realignment',
                          'maxrunning': 1000,
                          }

exe_reannotation = "/data/luojaa/proteinfer/proteinfer.py"
proteinfer_swarm_opts = {'threads-per-process': 8,
                          'gb-per-process': '10',
                          'time': '72:00:00',
                          'logdir': '/data/luojaa/log/swarm_out/',
                          'job-name': 'proteinfer',
                          'maxrunning': 1000,
                          }
clean_swarm_opts = {'threads-per-process': 8,
                          'gb-per-process': '25',
                          'time': '72:00:00',
                          'logdir': '/data/luojaa/log/swarm_out/',
                          'job-name': 'clean_lscratch',
                          'maxrunning': 1000,
                          'gres': 'lscratch:10',
                          }
# -

# hhsuite database formatting
hhconsensus_opts = '-M 50 -maxres 65535 -v 0'
# hhmake requires at least 2 shown sequences or else it segfaults...?
hhmake_opts = '-seq 2 -v 0'
cstranslate_opts = '-f -x 0.3 -c 4 -I a3m'

# hhsuite searches
exe_hhsuite_search = path_root + 'hhsuite_search.py'
#note that the alignment need to be present for parsing of full names using parse_HHsuite
hhblits_search_opts = '-n 1 -p 30 -z 1 -Z 100 -B 100 -b 1'
hhblits_swarm_opts = {'threads-per-process': 24,
                      'gb-per-process': '100',
                      'time': '48:00:00',
                      'gres': 'lscratch:100',
                      'logdir': '/data/luojaa/log/swarm_out/',
                      'job-name': 'hhsuite_search',
                      'maxrunning': 500,
                      'module': 'mmseqs,clustalo,muscle,python/3.9,hhsuite,mafft/7.475,IQTREE,FastTree'
                      }
# hhsuite merging/parsing search results
exe_hhsuite_parse = path_root + 'hhsuite_parse.py'
parse_swarm_opts = {'threads-per-process': 8,
                      'gb-per-process': '80',
                      'time': '8:00:00',
                      'logdir': '/data/luojaa/log/swarm_out/',
                      'job-name': 'hhsuite_parse',
                      'maxrunning': 50,
                      }

#more lscratch memory to accomodate uniref
hhblits_swarm_opts_uniref = {'threads-per-process': 24,
                      'gb-per-process': '30',
                      'time': '48:00:00',
                      'gres': 'lscratch:300',
                      'logdir': '/data/luojaa/log/swarm_out/',
                      'job-name': 'hhsuite_search',
                      'maxrunning': 500,
                      }

#iqtree saturates at 8 threads for most alignments
microcosm_format_opts = {
                         # 'original_query_DB': '/data/luojaa/eukgen/mmseqs/kog_proteins',
                         # 'original_target_DB': '/data/luojaa/eukgen/mmseqs/kog_proteins',
                         'original_query_DB': '/data/luojaa/eukgen/mmseqs_victor/euk72_ep/euk72_ep.repseq',
                         'original_target_DB': '/data/luojaa/eukgen/mmseqs_victor/prok2111_as/prok2111_as.repseq.minHGT',
                         # 'original_query_DB': '/data/luojaa/eukgen/processing/asgard2023/asgard2023_clean.repseq',
                         # 'original_target_DB': '/data/luojaa/eukgen/processing/prok2111/prok2111.repseq',
                         'taxonomy_mapping': merged_protein_tree_taxonomy,
                         'max_euk_sequences': 30,
                         'max_prok_sequences': 70,
                         'filter_entropy': 0.15,
                         'filter_length_frac': 0.2,
                         'threads': 8,
                         'muscle_reps': 3,
                         'muscle_timeout': 3600,
                         'evo_model_params': '-m MFP -mset LG,Q.pfam --cmin 4 --cmax 12'
                         }

exe_icarus_alignment = path_root + 'icarus_pairwise_align.py'
icarus_swarm_opts = {'threads-per-process': 16,
                     'gb-per-process': '75',
                     'time': '24:00:00',
                     'logdir': '/data/luojaa/log/swarm_out/',
                     'job-name': 'icarus_alignment',
                     'maxrunning': 1000}
