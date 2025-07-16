# python global packages
import pandas as pd
import subprocess
import os
from multiprocessing import current_process


# prepare microcosm file structure
def structure_files(root, query, query_hits, query_clusters, target_clusters):
    root_query = root+query+'/'
    os.mkdir(root_query)

    print(f'Q:{query}')
    query_members = query_clusters.loc[query, 'acc']
    hits = query_hits.loc[query]
    members = target_clusters.loc[hits]

    # write query accessions to .acc
    with open(root_query+f'{query}.query.acc', 'w') as outfile:
        outfile.write(pd.DataFrame(query_members).to_csv(sep='\t', header=None, index=None))

    with open(root_query+f'{query}.target.acc', 'w') as outfile:
        outfile.write(pd.DataFrame(members).to_csv(sep='\t', header=None))

# main script for formatting query and target .fasta and .cluster.tsv files

def prepare_mmseqs(file_root, file_basename, original_query_DB, original_target_DB, save_intermediate_files=False):
    import subprocess

    # configure paths and flags
    thread = 1
    basename = file_root + file_basename

    threadID_string = f'{thread} | {file_basename}:'
    # subprocess.run(f'mkdir {query_root}/tmp'.split())

    query_acc = f"{basename}.query.acc"
    target_acc = f"{basename}.target.acc"

    query_seqDB = f"{basename}.query.DB"
    target_seqDB = f"{basename}.target.DB"

    query_fasta = f"{basename}.query.fasta"
    target_fasta = f"{basename}.target.fasta"

    print(f'{threadID_string} Started \n', end='')

    # create a new DB for the query and target sequences
    print(f'{threadID_string} Preparing mmseqs data for query hits\n', end='')
    subprocess.run(f"mmseqs createsubdb -v 0 --id-mode 1 --subdb-mode 1 {query_acc} {original_query_DB} {query_seqDB}".split())
    subprocess.run(f'mmseqs convert2fasta -v 0 {query_seqDB} {query_fasta}'.split())

    print(f'{threadID_string} Preparing mmseqs data for target hits\n', end='')
    subprocess.run(f"mmseqs createsubdb -v 0 --id-mode 1 --subdb-mode 1 {target_acc} {original_target_DB} {target_seqDB}".split())
    subprocess.run(f'mmseqs convert2fasta -v 0 {target_seqDB} {target_fasta}'.split())

    # create tsv with lineage information from spoofed clustering result
    # create a file with double mmseqs_id column for DB to mimic result.tsv
    print(f'{threadID_string} Retreiving taxonomy lineages\n', end='')

    subprocess.run(f"cut -f1 {query_seqDB}.index | awk -v OFS='\t' '{{print $1, $1}}' > {query_seqDB}.tmp", shell=True)
    subprocess.run(f'mmseqs tsv2db -v 0 {query_seqDB}.tmp {query_seqDB}.tmp'.split())
    subprocess.run(f'mmseqs addtaxonomy -v 0 --tax-lineage 1 {query_seqDB} {query_seqDB}.tmp {query_seqDB}.tax'.split())
    subprocess.run(f'mmseqs createtsv -v 0 {query_seqDB} {query_seqDB} {query_seqDB}.tax {query_seqDB}.tsv'.split())

    subprocess.run(f"cut -f1 {target_seqDB}.index | awk -v OFS='\t' '{{print $1, $1}}' > {target_seqDB}.tmp",
                   shell=True)
    subprocess.run(f'mmseqs tsv2db -v 0 {target_seqDB}.tmp {target_seqDB}.tmp'.split())
    subprocess.run(
        f'mmseqs addtaxonomy -v 0 --tax-lineage 1 {target_seqDB} {target_seqDB}.tmp {target_seqDB}.tax'.split())
    subprocess.run(f'mmseqs createtsv -v 0 {target_seqDB} {target_seqDB} {target_seqDB}.tax {target_seqDB}.tsv'.split())

    # clean
    if not save_intermediate_files:
        subprocess.run(f'mmseqs rmdb -v 0 {query_seqDB}'.split())
        subprocess.run(f'mmseqs rmdb -v 0 {target_seqDB}'.split())
        subprocess.run(f'mmseqs rmdb -v 0 {query_seqDB}_h'.split())
        subprocess.run(f'mmseqs rmdb -v 0 {target_seqDB}_h'.split())
        subprocess.run(f'rm {target_seqDB}.tmp* {target_seqDB}.tax*', shell=True)
        subprocess.run(f'rm {query_seqDB}.tmp* {query_seqDB}.tax*', shell=True)
        subprocess.run(f"find {file_root} -type l -exec unlink {{}} \\;", shell=True)

    return query_fasta, target_fasta


# read fasta, align with FAMSA, filter and construct FastTree,
# crop eaves to size and write new fasta with only cropeed tree leaf sequences
# if supplied with taxonomy, will do taxonomically aware cropping
def fasta_reduce_size(base_fasta, threads, max_leaf_size, filter_entropy, 
                      save_intermediate_files=True, taxDF=None, min_clade_purity=0.9, min_clade_size=2,
                      predelete_outliers = False, leaf_cutoff = 0.95,  clade_cutoff = 0.95 
                     ):

    thread = current_process().pid
    threadID_string = f'{thread} | {base_fasta}:'

    from ete3 import Tree
    from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy
    from core_functions.tree_functions import get_outlier_nodes_by_lognorm, remove_outlier_nodes
    from paths_and_parameters import exe_fasttree, exe_famsa


    with open(base_fasta, 'r') as fastafile:
        size = fastafile.read().count('>')

    # threads for FastTree
    os.environ['OMP_NUM_THREADS'] = str(threads)

    famsa_logfile = open(f'{base_fasta}.famsa.log', 'a')
    fasttree_logfile = open(f'{base_fasta}.fasttree.log', 'a')

    # align seqs
    print(threadID_string + ' Aligning FAMSA')

    famsa_command = exe_famsa + f' -t {threads} {base_fasta} {base_fasta}.famsa'
    subprocess.run(famsa_command.split(), stdout=famsa_logfile, stderr=famsa_logfile)

    # backup original alignment
    if save_intermediate_files:
        subprocess.run(f'cp {base_fasta}.famsa {base_fasta}.famsa.b'.split())

    # filter euk by entropy
    print(threadID_string + ' Filtering by entropy')

    aln = fasta_to_dict(file=f'{base_fasta}.famsa')
    aln_filter = filter_by_entropy(aln, filter_entropy)
    dict_to_fasta(aln_filter, write_file=f'{base_fasta}.famsa')

    # construct a FastTree
    print(threadID_string + ' Constructing FastTree')

    fasttree_command = exe_fasttree + f" -gamma -out {base_fasta}.fasttree {base_fasta}.famsa"
    subprocess.run(fasttree_command.split(), stdout=fasttree_logfile, stderr=fasttree_logfile)
    
    tree = Tree(f'{base_fasta}.fasttree')
    
    # remove outlier leaves and nodes to avoid cropping to divergent set
    if predelete_outliers:
        
        print(f'{threadID_string} removing outliers before cropping')

        cropped_tree = tree.copy()
        outlier_nodes, dist_series, fits, cutoff = get_outlier_nodes_by_lognorm(cropped_tree, p_low=0, p_high=leaf_cutoff, only_leaves=True, deletion_cutoff=0.2)
        cropped_tree = remove_outlier_nodes(cropped_tree, outlier_nodes)
        print(f'After removing outliers leaves Tree has {len(cropped_tree)} leaves')

        outlier_nodes, dist_series, fits, cutoff = get_outlier_nodes_by_lognorm(cropped_tree, p_low=0, p_high=clade_cutoff, drop_leaves=True, only_leaves=False, deletion_cutoff=0.15)
        cropped_tree = remove_outlier_nodes(cropped_tree, outlier_nodes)
        print(f'After removing outliers Tree has {len(cropped_tree)} leaves')
        
        tree = cropped_tree.copy()

    # reduce leaves to size
    print(threadID_string + ' Cropping Leaves')

        
    # if taxonomic information is supplied crop while optimizing taxa
    # retain at least minimum clade size for monophyletic taxa
    if taxDF is not None:
        from core_functions.tree_functions import crop_leaves_to_size_considering_taxa
        tree, leafDF, crop_dict = crop_leaves_to_size_considering_taxa(tree, taxDF, max_leaf_size,
                                                                       min_clade_size=min_clade_size, 
                                                                       min_clade_purity=min_clade_purity,
                                                                       LCA_search_depth=5)

        # save leaf mappings from DF
        print(threadID_string + f' Writing .leaf_mapping after cropping')
        leafDF.to_csv(f'{base_fasta}.leaf_mapping', header=None, sep='\t')

    # if no taxonomy, crop leaves as per normal
    else:
        from core_functions.tree_functions import crop_leaves_to_size
        tree, crop_dict = crop_leaves_to_size(tree, max_leaf_size)

        #save cropped leaf mappings from dict
        with open(base_fasta + '.leaf_mapping', 'w') as outfile:
            for leaf, cropped in crop_dict.items():
                for acc in cropped:
                    outfile.write(f'{acc}\t{leaf}\n')

    # write cropped fasta
    cropped_leaves = tree.get_leaf_names()
    cropped_aln = {key: value.replace('-', '') for key, value in aln.items() if key in cropped_leaves}

    # replace original fasta with cropped version and save backup
    if save_intermediate_files:
        subprocess.run(f'cp {base_fasta} {base_fasta}.uncropped'.split())
    dict_to_fasta(cropped_aln, write_file=base_fasta)

    famsa_logfile.close()
    fasttree_logfile.close()

    if not save_intermediate_files:
        delete_files = f'rm {base_fasta}.famsa.log {base_fasta}.fasttree.log {base_fasta}.famsa {base_fasta}.fasttree {base_fasta}.leaf_mapping'
        subprocess.run(delete_files.split())

    return base_fasta


# main analysis function for tree annotation
# assigns taxonomy to leaves, identifies and trims outliers, assigns LCA nodes and performs distance calculations
# outputs a treeDF with tabulated data for each eukaryotic clade and its corresponding prokaryotic sister-clades
def tree_analysis(tree_file, leaf_mapping, tree_name,
                  outlier_CDF_low=0,
                  outlier_CDF_high=0.99,
                  tree_crop_cutoff=0.30,
                  delete_outliers=True,
                  prok_clade_size=3,
                  prok_clade_purity=0.80,
                  euk_clade_size=5,
                  euk_clade_purity=0.8,
                  exclude_nested_LCAs=True,
                  consider_closest_n_prok_LCAs=12):

    import pandas as pd
    from ete3 import Tree
    from core_functions.tree_functions import get_outlier_nodes_by_lognorm, map_leafDF, \
        get_multiple_soft_LCA_by_relative_purity, weighted_midpoint_root

    # load tree
    print(f'Reading Tree from {tree_file}')
    tree = Tree(tree_file)

    # load leaf mapping from tree cropping
    print(f'Loading leaf data from {leaf_mapping}')
    leafDF = pd.read_csv(leaf_mapping, sep='\t', names=['acc', 'leaf', 'superkingdom', 'full_class', 'rank'],
                         index_col=0)

    # merge all euk classification in to 'Eukaryota' save old class
    leafDF['class'] = leafDF.full_class
    leafDF.loc[leafDF[leafDF.superkingdom == 'Eukaryota'].index, 'class'] = 'Eukaryota'

    # reroot at midpoint before outlier calculation
    tree = weighted_midpoint_root(tree)

    # outlier detection
    print(f'Running outlier detection')
    outlier_nodes, dist_series, fit_params, cutoffs = get_outlier_nodes_by_lognorm(tree, outlier_CDF_low,
                                                                                   outlier_CDF_high,
                                                                                   deletion_cutoff=tree_crop_cutoff)

    # update leafDF with deleted nodes
    outlier_leaves = [leaf.name for node in outlier_nodes for leaf in node.get_leaves()]
    leafDF.loc[outlier_leaves, 'leaf'] = 'DELETED'

    # detach all outlier clades
    if delete_outliers:
        for node in outlier_nodes:
            node.detach()

    # repair tree and collapse unifurcations
    # accoring to the API this shouldn't be neccesary but the detach() funtionality appears inconsistent
    for node in tree.traverse():
        if len(node.children) == 1 or (node.is_leaf() and node.name == ''):
            node.delete()

    # annotate tree with taxonomic info accounting for collapsed branches
    tree = map_leafDF(tree, leafDF)

    # find the best bacterial soft LCA for all taxa
    taxa_list = leafDF[(leafDF.superkingdom.isin(['Bacteria', 'Archaea'])) & (leafDF.leaf != 'DELETED')][
        'class'].unique()

    lca_dict = {}
    for taxa in sorted(taxa_list):
        lca_dict[taxa] = get_multiple_soft_LCA_by_relative_purity(tree, taxa,
                                                                  n_best=9999,
                                                                  min_size=prok_clade_size,
                                                                  min_purity=prok_clade_purity,
                                                                  max_depth=10)

    # find all euk sof LCAs meeting criteria
    lca_dict['Eukaryota'] = get_multiple_soft_LCA_by_relative_purity(tree, 'Eukaryota',
                                                                     n_best=9999,
                                                                     min_size=euk_clade_size,
                                                                     min_purity=euk_clade_purity,
                                                                     max_depth=10)

    # display LCA results
    accepted_taxa = [taxa for taxa, node in lca_dict.items() if node != [] and taxa != 'Eukaryota']
    rejected_taxa = [taxa for taxa, node in lca_dict.items() if node == [] and taxa != 'Eukaryota']

    if accepted_taxa:
        print(f'Found acceptable LCAs with size > {prok_clade_size} and purity {prok_clade_purity} for: ')
        for taxa in accepted_taxa:
            node_data = lca_dict[taxa]
            for node in sorted(node_data, key=lambda x: -x[1]):
                print(f'    {taxa + ":":<25} size: {node[1]}\t weight: {round(node[3], 2)}')
                pass

        print(f'\nNo acceptable LCAs with size > {prok_clade_size} and purity {prok_clade_purity} for:')
        print(*sorted(rejected_taxa), sep='\n')

    else:
        print(f'WARNING! Found no acceptable LCAs for Prokayotes, will exit.')
        pass

    print(f'\nFound {len(lca_dict["Eukaryota"])} LCA nodes for Eukaryota')
    for node_data in lca_dict['Eukaryota']:
        print(f'    {"Eukarya:":<25} size: {node_data[1]}\t weight: {round(node_data[3], 2)}')
        pass
    print()

    # if more than three good euk nodes are identified warn and restruct analysis
    if len(lca_dict['Eukaryota']) > 3:
        print(
            f'WARNING! High paraphyly in Eukaryota, considering the three largest clades of {len(lca_dict["Eukaryota"])} total!')
        lca_dict['Eukaryota'] = lca_dict['Eukaryota'][:3]

    # revert to original eukaryotic classes on tree
    leafDF.columns = ['leaf', 'superkingdom', 'class', 'rank', 'filter_class']
    tree = map_leafDF(tree, leafDF)

    # calculate distances and format DataFrame
    euk_dist_dict = {}

    # separate euk and prok LCA nodes
    all_LCA_data = [node for taxa, nodes in lca_dict.items() for node in nodes]
    all_LCA_nodes = [node[0] for node in all_LCA_data]
    euk_LCA_nodes = [node for node in all_LCA_data if node[0].lca_taxa == 'Eukaryota']
    prok_LCA_nodes = [node for node in all_LCA_data if node[0].lca_taxa != 'Eukaryota']

    # for all eukaryotic nodes
    for i, euk_node_data in enumerate(euk_LCA_nodes):

        # store node reference and basic information
        euk_node = euk_node_data[0]
        euk_node_size = euk_node_data[1]
        euk_clade_is_leaf = euk_node.is_leaf()
        euk_node_acc = euk_node.get_closest_leaf()[0].name

        # get the LCA for the identified EUK node as well as the unique class members
        euk_node_LCA, euk_node_scope, euk_node_scope_len = retrieve_LCA_scope(euk_node.get_leaf_names(), leafDF)

        # for all the best prok nodes per taxa
        for j, prok_node_data in enumerate(prok_LCA_nodes):

            # store node reference and basic information
            prok_node = prok_node_data[0]
            prok_node_acc = prok_node.get_closest_leaf()[0].name
            prok_node_is_clade = prok_node.is_leaf()

            # calculate distance metrics
            # absolute node-to-node distance
            dist = euk_node.get_distance(prok_node)
            top_dist = euk_node.get_distance(prok_node, topology_only=True)

            # gabaldon stem length 2016 calculation
            lca_node_prok_euk = euk_node.get_common_ancestor(prok_node)
            raw_stem_length = lca_node_prok_euk.get_distance(euk_node)

            # cannot normalize if clade is leaf
            if euk_clade_is_leaf:
                stem_length = raw_stem_length
                median_euk_branch_length = -1

            # otherwise normalise as per gabaldon
            else:
                euk_branch_lengths = pd.Series([euk_node.get_distance(leaf) for leaf in euk_node.get_leaves()])
                median_euk_branch_length = euk_branch_lengths.median()
                stem_length = raw_stem_length / median_euk_branch_length

            # save data as temp_index: treename_"clade_number", clade_name, clade_size, prok_clade_name, prok_taxa...
            # final value of False is for later filtering purposes
            euk_dist_dict[str(i) + '_' + str(j)] = [tree_name, euk_node_acc, euk_node_data[1], euk_node_data[3],
                                                    euk_clade_is_leaf, euk_node_LCA, euk_node_scope, euk_node_scope_len,
                                                    prok_node_acc, prok_node_data[1], prok_node_data[3],
                                                    prok_node_is_clade, prok_node.lca_taxa, dist, top_dist,
                                                    raw_stem_length, median_euk_branch_length, stem_length,
                                                    euk_node, prok_node, False]

    tree_data = pd.DataFrame.from_dict(euk_dist_dict, orient='index',
                                       columns=['tree_name', 'euk_clade_rep', 'euk_clade_size', 'euk_clade_weight',
                                                'euk_leaf_clade', 'euk_LCA', 'euk_scope', 'euk_scope_len',
                                                'prok_clade_rep', 'prok_clade_size', 'prok_clade_weight',
                                                'prok_leaf_clade', 'prok_taxa', 'dist', 'top_dist',
                                                'raw_stem_length', 'median_euk_leaf_dist', 'stem_length',
                                                'euk_node_ref', 'prok_node_ref', 'include'])

    # if no euk clades were found return empty dataframe
    if tree_data.shape[0] == 0:
        return tree, tree_data

    # drop temporary index
    tree_data.sort_values(by=['top_dist', 'dist'], ascending=True, inplace=True)
    tree_data.set_index('tree_name', inplace=True)

    # check if LCA is decendants of any LCA other node
    # this CURRENTLY breaks downstream formatting of contraints analysis
    if exclude_nested_LCAs:
        print(f'Filtering to remove nested LCAs, nodes that are children of other nodes')
        nested_nodes = []
        for test_node in all_LCA_nodes:

            if any(node in all_LCA_nodes for node in test_node.get_ancestors()):
                test_node_acc = test_node.get_closest_leaf()[0].name
                # print(f'EXCLUDING LCA node {test_node_acc} for {test_node.lca_taxa}')
                nested_nodes.append(test_node)

        # exclude nested LCAs
        tree_data['decendant'] = [any(node in nested_nodes for node in data) for data in
                                  tree_data[['euk_node_ref', 'prok_node_ref']].values]
        tree_data = tree_data[~(tree_data.decendant)]

        print()

    # loop through data including only n closest prok LCAs per euk LCAs
    print(
        f'Consider only closest clades for downstream analysis, clade limit per Eukaryote clade is {consider_closest_n_prok_LCAs}')

    for euk_clade_rep, data in tree_data.groupby('euk_clade_rep'):
        print(f'Including the following LCA list for Eukaryotic LCA: {euk_clade_rep}')
        data.iloc[0].euk_node_ref.add_feature('LCA', 'Eukaryota')

        # consider first n prok LCAs by slicing sorted list
        for i, prok_data in data.iloc[:consider_closest_n_prok_LCAs].iterrows():
            print(
                f'   {prok_data.prok_taxa + ":":<30} {prok_data.prok_clade_rep:<25} top_dist: {prok_data.top_dist:<8}')
            prok_data.prok_node_ref.add_feature('LCA', prok_data.prok_taxa)

            # add a reference to keep all prok LCAs, quite inefficient
            tree_data.loc[(tree_data.euk_clade_rep == euk_clade_rep) &
                          (tree_data.prok_clade_rep == prok_data.prok_clade_rep), 'include'] = True

        print()

    # delete superflous LCA data and node references
    tree_data = tree_data[tree_data.include]
    tree_data.drop(['euk_node_ref', 'prok_node_ref', 'include', 'decendant'], axis=1, inplace=True)

    tree.ladderize()

    return tree, tree_data


# get LCA ancestor for post tree partitioned EUK clade
def retrieve_LCA_scope(LCA_leaf_names, leafDF):
    # hard coded data from ete3 for faster retrevial
    # modified to more finely sample eukarya, additions of "Diaphoretickes" and "Amphorea", "TSAR", "Excavata" etc.
    # for clades with orphans like "opistokonta" which are superclades of other clades like "metazoa" I duble the last entry to avoid popping from an empty list
    class_lineages = {'Cyanobacteriota': ['root', 'cellular organisms', 'Bacteria', 'Terrabacteria group', 'Cyanobacteriota/Melainabacteria group', 'Cyanobacteriota'],
                      'Gammaproteobacteria': ['root', 'cellular organisms', 'Bacteria', 'Pseudomonadota', 'Gammaproteobacteria'],
                      'Bacillota': ['root', 'cellular organisms', 'Bacteria', 'Terrabacteria group', 'Bacillota'],
                      'Rhodophyta': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'Archaeplastida','Rhodophyta'],
                      'Dinophyceae': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar','Alveolata', 'Dinophyceae'],
                      'Cryptophyceae': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes','Cryptista' 'Cryptophyceae'],
                      'Euglenida': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Excavata', 'Discoba','Euglenozoa', 'Euglenida'],
                      'Chlorophyta': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'Archaeplastida','Viridiplantae', 'Chlorophyta'],
                      'Oomycota': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar','Stramenopiles', 'Oomycota'],
                      'Ascomycota': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Fungi', 'Dikarya', 'Ascomycota'],
                      'Basidiomycota': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Fungi', 'Dikarya', 'Basidiomycota'],
                      'Heterolobosea': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Excavata', 'Discoba', 'Heterolobosea'],
                      'Apicomplexa': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar', 'Alveolata', 'Apicomplexa'],
                      'Ciliophora': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar','Alveolata', 'Ciliophora'],
                      'Microsporidia': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Fungi', 'Fungi incertae sedis', 'Microsporidia'],
                      'Choanoflagellata': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta', 'Choanoflagellata'],
                      'Alphaproteobacteria': ['root', 'cellular organisms', 'Bacteria', 'Pseudomonadota','Alphaproteobacteria'],
                      'Betaproteobacteria': ['root', 'cellular organisms', 'Bacteria', 'Pseudomonadota', 'Betaproteobacteria'],
                      'Deltaproteobacteria': ['root', 'cellular organisms', 'Bacteria', 'delta/epsilon subdivisions', 'Deltaproteobacteria'],
                      'Euryarchaeota': ['root', 'cellular organisms', 'Archaea', 'Euryarchaeota'],
                      'Asgard': ['root', 'cellular organisms', 'Archaea', 'Asgard'],
                      'Campylobacterota': ['root', 'cellular organisms', 'Bacteria', 'Campylobacterota'],
                      'Opisthokonta': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Opisthokonta'],
                      'Metazoa': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Metazoa'],
                      'Streptophyta': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'Archaeplastida','Viridiplantae', 'Streptophyta'],
                      'Glaucocystophyceae': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes','Archaeplastida', 'Glaucocystophyceae'],
                      'Acidobacteriota': ['root', 'cellular organisms', 'Bacteria', 'Acidobacteriota'],
                      'Ancyromonadidae': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Ancyromonadida', 'Ancyromonadidae'],
                      'Fungi incertae sedis': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa','Opisthokonta', 'Fungi', 'Fungi incertae sedis', 'Fungi incertae sedis'],
                      'Ichthyosporea': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Ichthyosporea'],
                      'Malawimonadidae': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Excavata','Malawimonadida', 'Malawimonadidae'],
                      'Thermoplasmata': ['root', 'cellular organisms', 'Archaea', 'Candidatus Thermoplasmatota','Thermoplasmata'],
                      'Thermotogae': ['root', 'cellular organisms', 'Bacteria', 'Thermotogota', 'Thermotogae'],
                      'Chloroflexota': ['root', 'cellular organisms', 'Bacteria', 'Terrabacteria group','Chloroflexota'],
                      'Thermodesulfobacteriota': ['root', 'cellular organisms', 'Bacteria', 'Thermodesulfobacteriota'],
                      'Actinomycetota': ['root', 'cellular organisms', 'Bacteria', 'Terrabacteria group','Actinomycetota'],
                      'Chlamydiia': ['root', 'cellular organisms', 'Bacteria', 'PVC group', 'Chlamydiota','Chlamydiia'],
                      'Rhizaria': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar', 'Rhizaria'],
                      'Mycoplasmatota': ['root', 'cellular organisms', 'Bacteria', 'Terrabacteria group','Mycoplasmatota'],
                      'Apusozoa': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Apusozoa'],
                      'Amoebozoa': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Amoebozoa'],
                      'Jakobida': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Excavata', 'Discoba', 'Jakobida'],
                      'Telonemida': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Telonemida'],
                      'Colpodellida': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar','Alveolata', 'Colpodellida'],
                      'Tsukubamonadidae': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Excavata', 'Discoba','Tsukubamonadida', 'Tsukubamonadidae'],
                      'Breviatea': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Breviatea'],
                      'Colponemidia': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar', 'Alveolata', 'Colponemida', 'Colponemidia'],
                      'Aphelidea': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta', 'Aphelida', 'Aphelidea'],
                      'PVC group': ['root', 'cellular organisms', 'Bacteria', 'PVC group'],
                      'FCB group': ['root', 'cellular organisms', 'Bacteria', 'FCB group'],
                      'TACK group': ['root', 'cellular organisms', 'Archaea', 'TACK group'],
                      'DPANN group': ['root', 'cellular organisms', 'Archaea', 'DPANN group'],
                      'Fonticulaceae': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta','Rotosphaerida', 'Fonticulaceae'],
                      'Hemimastigophora': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'Hemimastigophora'],
                      'Perkinsozoa': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar','Alveolata', 'Perkinsozoa'],
                      'Rhodelphea': ['root', 'cellular organisms', 'Eukaryota', 'Rhodelphea'],
                      'Haptista': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'Haptista'],
                      'CRuMs': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'CRuMs'],
                      'Metamonada': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Excavata', 'Metamonada'],
                      'Eukaryota incertae sedis': ['root', 'cellular organisms', 'Eukaryota','Eukaryota incertae sedis'],
                      'Bigyra': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar', 'Stramenopiles', 'Bigyra'],
                      'Filasterea': ['root', 'cellular organisms', 'Eukaryota', 'Amorphea', 'Obazoa', 'Opisthokonta', 'Filasterea'],
                      'Ochrophyta': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar', 'Stramenopiles', 'Ochrophyta'],
                      'Sar': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'TSar', 'Sar', 'Sar'],
                      'Prasinodermophyceae': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes','Archaeplastida', 'Viridiplantae', 'Prasinodermophyta', 'Prasinodermophyceae'],
                      'Myxococcota': ['root', 'cellular organisms', 'Bacteria', 'Myxococcota'],
                      'Nebulidea': ['root', 'cellular organisms', 'Eukaryota', 'Diaphoretickes', 'Provora', 'Nebulidia','Nebulidea'],
                      'Bdellovibrionota': ['root', 'cellular organisms', 'Bacteria', 'Bdellovibrionota'],
                      'Spirochaetota': ['root', 'cellular organisms', 'Bacteria', 'Spirochaetota'],
                      'Deinococcota': ['root', 'cellular organisms', 'Bacteria', 'Terrabacteria group', 'Deinococcota'],
                      'Nitrospirota': ['root', 'cellular organisms', 'Bacteria', 'Nitrospirota'],
                      'DELETE': ['DELETE']}

    leaf_clades = leafDF[leafDF.index.isin(LCA_leaf_names) & (leafDF['class'] != 'DELETE')]['class'].unique()

    # if no class names returned something is wrong
    if len(leaf_clades) < 1:
        return 'ERROR', 'ERROR', 0

    # if the group is monophyletic return early as LCA and scope is itself
    if len(leaf_clades) == 1:
        LCA_clade = leaf_clades[0]

        return LCA_clade, LCA_clade, 1

    lineages = [class_lineages[clade].copy() for clade in leaf_clades]

    print(lineages)

    # pop first element of all lists until set of the elements is not single
    LCA_clade = ''
    while True:
        top_clade = {lin.pop(0) for lin in lineages}
        if len(top_clade) == 1:
            LCA_clade = top_clade.pop()

        else:
            break

    print(LCA_clade)

    return LCA_clade, '|'.join(leaf_clades), len(leaf_clades)


# generate a set of constraint trees for each euk_LCA * prok_LCA pair for IQtree2 -g constraint testing
# ((euk_LCA, prok_sister),(remaining_prok));
# return a dataframe of jobs to be run
# clade size dictates the maximum number of sequences forming the contrained clades
# CURRENLTY ONLY ALLOWING FOR VARIABLE EUK CLADE SIZE ONLY USING ONE PROK SEQUENCE
def format_constraint_analysis(root, basename, tree_data, clade_size=10):

    from core_functions.helper_functions import fasta_to_dict, dict_to_fasta
    import pandas as pd
    from ete3 import Tree

    # initialise filepaths
    treefile = root + basename + '.merged.fasta.muscle.treefile.annot'
    msafile = root + basename + '.merged.fasta.muscle'
    constraint_base = root + '/constraint_analysis/'

    # create analysis output folder
    subprocess.run(f'mkdir {constraint_base}'.split())

    # find all annotated LCA nodes
    tree = Tree(treefile)
    euk_LCA_nodes = []
    prok_LCA_nodes = []

    for node in tree.traverse():
        if 'LCA' in node.features:
            if node.LCA == 'Eukaryota':
                euk_LCA_nodes.append(node)
            else:
                prok_LCA_nodes.append(node)


    # extract all contrstaint tree pairs from tree_data as not all EUK LCAs are to be compared to all prok LCAs
    constraint_pairs = tree_data[['euk_clade_rep', 'prok_clade_rep']].values.tolist()

    # initialize values for jobDF creation
    constraint_data = {}
    index = 0

    # evaluate all euk_LCA nodes individually against all annotated prok_LCAs
    for euk_node in euk_LCA_nodes:

        # node representative is closest_leaf in clade
        euk_node_acc = euk_node.get_closest_leaf()[0].name
        euk_accs = {node for node in euk_node.get_leaf_names()[:clade_size]}

        # output file for filtered alignment saving
        msafile_out = f'{constraint_base}/{euk_node_acc}_constraint.fasta.muscle'

        for n, prok_node in enumerate(prok_LCA_nodes):
            # node representative is closest_leaf in clade
            prok_node_acc = prok_node.get_closest_leaf()[0].name

            # if current euk prok pair not in tree_data pairs exclude it
            if [euk_node_acc, prok_node_acc] not in constraint_pairs:
                print(f'Constraint "{euk_node_acc}, {prok_node_acc}" is excepted as pair does not exist in tree data')
                continue

            # sister group and outgroups
            prok_sister_accs = {name for name in prok_node.get_leaf_names()}

            # outgroup is all prok LCA sequences which are not in the siste rset or a decendant of the euk LCA
            additional_prok_accs = {name for node in prok_LCA_nodes for name in node.get_leaf_names() if
                                    name not in prok_sister_accs.union(euk_accs)}

            # handle cases of singleton outgroup branches
            # ((a,b,c),(d)) should be ((a,b,c),d)

            if len(additional_prok_accs) < 2:
                constraint = f"(({','.join(euk_accs)},{','.join(prok_sister_accs)}),{','.join(additional_prok_accs)});"

            else:
                constraint = f"(({','.join(euk_accs)},{','.join(prok_sister_accs)}),({','.join(additional_prok_accs)}));"

            constraint_file = f'{constraint_base}/{euk_node_acc}_{prok_node_acc}.constraint'
            with open(constraint_file, 'w') as tree:
                tree.write(constraint)

            # save relevant paths and tax info for downstream processing
            constraint_data[index] = [euk_node_acc]
            constraint_data[index].append(prok_node_acc)
            constraint_data[index].append(prok_node.LCA)
            constraint_data[index].append(constraint_file)
            constraint_data[index].append(msafile_out)
            index += 1

        # generate list of all sequences for MSA filtering
        all_prok_accs = {name for node in prok_LCA_nodes for name in node.get_leaf_names()}
        all_LCA_accs = list(euk_accs.union(all_prok_accs))

        # filter MSA for constraint analysis keeping only leaves in relevant LCAs
        aln = fasta_to_dict(file=msafile)
        filter_aln = {key: aln[key] for key in all_LCA_accs}
        dict_to_fasta(filter_aln, write_file=msafile_out, verbose=False)

    # format jobDF for strting IQtree runs
    column_names = ['euk_clade_rep', 'prok_clade_rep', 'prok_taxa', 'constraint_tree', 'constraint_msa']
    constraint_job_data = pd.DataFrame.from_dict(constraint_data, orient='index', columns=column_names)

    # nota all EUK LCAs are to be compared against all prok LCAs
    # filter contraint trees to match tree_data

    constraint_job_data.to_csv(f'{constraint_base}/constraint_trees.tsv', sep='\t', index=None)

    print(f'Formatted {constraint_job_data.shape[0]} guide trees for {len(euk_LCA_nodes)} eukaryotic LCAs')

    return constraint_job_data.set_index('euk_clade_rep')


# parse output from iqtree2 -z tree comparison result.iqtree file
# attempts to slice from "USER TREES" to "deltaL  :", rather fragile
# tree_names is list of names for each tree typically hypothesis test
def parse_IQtree_z_output(iqtreefile, tree_names=None):
    # read all lines
    with open(iqtreefile, 'r') as iqfile:
        raw_file = iqfile.read()

    # crop out table
    iqtree_data = raw_file.split('USER TREES\n')[1].split('deltaL  :')[0].split('\n')[6:-2]

    # format backup names if no tree_names
    if tree_names is None:
        tree_names = ['tree_' + str(i) for i in range(len(iqtree_data))]

    # split data lines and format into dict for pandas
    iqtree_dict = {}
    for i, line in enumerate(iqtree_data):
        data_list = line.split()
        iqtree_dict[tree_names[i]] = [i for i in data_list[1:]]

    # format pandasDF
    column_names = ['logL', 'deltaL', 'bp-RELL', 'bp-RELL_accept', 'p-KH', 'p-KH_accept', 'p-SH', 'p-SH_accept',
                    'c-ELW', 'c-ELW_accept', 'p-AU', 'p-AU_accept']
    test_data = pd.DataFrame.from_dict(iqtree_dict, orient='index', columns=column_names)

    # set numerical datatypes for relevant columns
    for col in ['logL', 'deltaL', 'bp-RELL', 'p-KH', 'p-SH', 'c-ELW', 'p-AU']:
        test_data[col] = pd.to_numeric(test_data[col])

    test_data = test_data.reset_index(names='constraint_tree_id')
    return test_data


# execute all constrained iqtree runs serially
def run_constraint_analysis(constraint_job_data, evo_model, threads):
    from paths_and_parameters import exe_iqtree

    # construct one constrained tree for each euk_LCA*prok_LCA pair
    for index, row in constraint_job_data.iterrows():
        alignment = row.constraint_msa
        constraint_tree = row.constraint_tree

        print(f'Running IQtree2 constrained tree analysis for {constraint_tree}')
        with open(f'{constraint_tree}.log', 'a') as iqtree_logfile:
            # run 1000 ultrafast bootstraps use given model, add -bnni for UFBoot model violations
            iqtree_command = f'{exe_iqtree} -s {alignment} -m {evo_model} -g {constraint_tree} --prefix {constraint_tree} --threads {threads} -B 1000 --redo'
            subprocess.run(iqtree_command.split(), stdout=iqtree_logfile, stderr=iqtree_logfile)

    # concatenate all resulting treefiles for tree testing per euk_LCA
    test_data_list = []

    for euk_LCA, job_data in constraint_job_data.groupby('euk_clade_rep'):
        alignment = job_data.constraint_msa.values[0]

        # maintain DF order for bash cat as iqtree -z follows input order
        bash_cat_order = ' '.join([treefile + '.treefile' for treefile in job_data.constraint_tree])
        forest_file = f'{alignment}.forestfile'
        subprocess.run(f"cat {bash_cat_order} > {forest_file}", shell=True)

        # run iqtree tree evaluation
        with open(f'{alignment}.log', 'a') as iqtree_logfile:
            print(f'Concatenating {job_data.shape[0]} trees for {euk_LCA} and evaluating')
            iqtree_command = f'{exe_iqtree} -s {alignment} -m {evo_model} -z {forest_file} -n 0 -zb 10000 -au --redo'
            subprocess.run(iqtree_command.split(), stdout=iqtree_logfile, stderr=iqtree_logfile)

        # format output database and append to data_list
        iqtreefile = alignment + '.iqtree'
        test_data = parse_IQtree_z_output(iqtreefile, tree_names=None)
        test_data['euk_clade_rep'] = [euk_LCA] * test_data.shape[0]
        test_data['prok_clade_rep'] = job_data.prok_clade_rep.values
        test_data['prok_taxa'] = [taxa for taxa in job_data['prok_taxa']]
        test_data_list.append(test_data)

    # merge all data lists and return
    all_test_data = pd.concat(test_data_list)
    all_test_data.drop('constraint_tree_id', axis=1, inplace=True)

    return all_test_data

# applies clade markings and styles to tree from tree_analysis:
def color_tree(tree, savefile=None, view_in_notebook=False):
    from ete3 import TreeStyle, NodeStyle, TextFace

    # fixes some rendering crasher as per
    # https://github.com/etetoolkit/ete/issues/296
    import os
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'

    # AVOID OVERWRITING OLD TREE
    annot_tree = tree.copy()

    ts = TreeStyle()
    # ts.title.add_face(TextFace(tree_header, fsize=8), column=0)
    ts.mode = 'r'
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.optimal_scale_level = ''
    ts.allow_face_overlap = True

    # set individual node styles
    default_node_style = NodeStyle()
    default_node_style['size'] = 0
    default_node_style['fgcolor'] = 'Black'

    LCA_euk_node_style = NodeStyle()
    LCA_euk_node_style['size'] = 10
    LCA_euk_node_style['fgcolor'] = '#5c7fe1'
    LCA_euk_node_style['bgcolor'] = '#eff3ff'
    
    LCA_asg_node_style = NodeStyle()
    LCA_asg_node_style['size'] = 10
    LCA_asg_node_style['fgcolor'] = '#dbcd6e'
    LCA_asg_node_style['bgcolor'] = '#fffeef'

    LCA_prok_node_style = NodeStyle()
    LCA_prok_node_style['size'] = 8
    LCA_prok_node_style['fgcolor'] = 'Black'
    LCA_prok_node_style['bgcolor'] = 'LightGray'

    for node in annot_tree.traverse():
        node.set_style(default_node_style)

        if node.is_leaf():
            node.add_face(TextFace(node.name, fsize=8), column=1)
            if 'taxa' in node.features:
                node.add_face(TextFace(' ' + str(node.taxa), fsize=8), column=2)
                node.add_face(TextFace(' ' + str(node.counts), fsize=8), column=3)

        if 'LCA' in node.features:
            if node.LCA == 'Eukaryota':
                node.set_style(LCA_euk_node_style)
                #node.add_face(TextFace(node.LCA, fsize=8), column = 1)
                
            elif node.LCA == 'Asgard':
                node.set_style(LCA_asg_node_style)

            else:
                node.set_style(LCA_prok_node_style)

    annot_tree.ladderize()

    if savefile != None:
        # write to pdf
        annot_tree.render(savefile + '.pdf', tree_style=ts)

    if view_in_notebook:

        tmp_file = f'./tmp/tmp_tree.png'
        _ = annot_tree.render(tmp_file, tree_style=ts)

        from IPython.display import Image
        return annot_tree, Image(tmp_file)

    return annot_tree


# try to merge two sets of fastas and construct a tree
# evaluate the grouping of original set labels the tree to see if the best partition is non-random
# runs famsa+FastTree cropping columns inbetween
def test_anneal_fastas(fastas, threads=8, filter_entropy=0, max_leaves=500, random_label_samples=200,
                      save_intermediate_files=False):
    
    from ete3 import Tree
    from Bio import SeqIO
    from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy
    from core_functions.tree_functions import crop_leaves_to_size, weighted_midpoint_root
    from paths_and_parameters import exe_fasttree, exe_famsa

    # threads for FastTree
    os.environ['OMP_NUM_THREADS'] = str(threads)

    merged_fasta = f'{fastas[0].rsplit("/", 1)[0]}/{"_vs_".join(fasta.split("/")[-1] for fasta in fastas)}'

    all_seqs = []
    leaf_data = []
    for fasta in fastas:
        name = fasta.split('/')[-1]
        seqs = [*SeqIO.parse(fasta, format='fasta')]
        leaf_data.append(pd.DataFrame({'acc':[s.id for s in seqs], 'set':[name]*len(seqs)}))
        all_seqs.extend(seqs)

    # recoed origin set of leaves and save to file
    leaf_mapping = pd.concat(leaf_data).set_index('acc')
    SeqIO.write(all_seqs, merged_fasta, 'fasta')

    print(f'Evaluating merger of {len(all_seqs)} sequences from {len(fastas)} sets')


    
    # align filter and construct tree for seqs
    famsa_logfile = open(f'{merged_fasta}.famsa.log', 'a')
    fasttree_logfile = open(f'{merged_fasta}.fasttree.log', 'a')

    famsa_command = exe_famsa + f' -t {threads} {merged_fasta} {merged_fasta}.famsa'
    subprocess.run(famsa_command.split(), stdout=famsa_logfile, stderr=famsa_logfile)

    aln = fasta_to_dict(file=f'{merged_fasta}.famsa')
    aln_filter = filter_by_entropy(aln, filter_entropy)
    dict_to_fasta(aln_filter, write_file=f'{merged_fasta}.famsa')

    fasttree_command = exe_fasttree + f" -gamma -out {merged_fasta}.fasttree {merged_fasta}.famsa"
    subprocess.run(fasttree_command.split(), stdout=fasttree_logfile, stderr=fasttree_logfile)

    
    
    # reduce leaves to size for faster bootstrapping
    tree = Tree(f'{merged_fasta}.fasttree')

    # annotate leaves with real set membership for calculation of real partition
    for i, l in enumerate(tree.get_leaves()):
        l.add_feature('set', leaf_mapping.loc[l.name,'set'][0])

    # reroot and save original tree
    tree = weighted_midpoint_root(tree)
    cropped_tree = tree.copy()

    # crop tree and resample leaf mapping
    print(f'Cropping joint tree to {max_leaves}')
    cropped_tree, crop_dict = crop_leaves_to_size(cropped_tree, max_size = max_leaves, 
                                          save_cropped_references=True, monophyletic=False, crop_dict=None)

    cropped_leaf_mapping = pd.DataFrame([[l.name, leaf_mapping.loc[l.name,'set']] for l in cropped_tree.get_leaves()], 
                                        columns = ['acc', 'set'])

    
    # bootstrap samples for evaluating partitioning score
    data = []
    for n in range(random_label_samples):

        # reassign random leaf mapping from cropped labels to cropped tree leaves
        rand_leaf_mapping = cropped_leaf_mapping.sample(frac=1, replace=False)['set'].values
        for i, l in enumerate(cropped_tree.get_leaves()):
                l.add_feature('set', rand_leaf_mapping[i])

        # find score of best partition over the randomized tree labels
        node = find_optimal_tree_partition(cropped_tree, attribute='set')
        data.append(node.score)
        
    data = pd.Series(data)
    dist = stats.expon(*stats.expon.fit(data))
    
    
    
    # annotate leaves with real set membership for calculation of real partition
    for i, l in enumerate(cropped_tree.get_leaves()):
        l.add_feature('set', leaf_mapping.loc[l.name,'set'])

    node = find_optimal_tree_partition(cropped_tree, attribute='set')

    p = (data[data >= node.score].count())/random_label_samples
    print(f'Random score outperformed real score ({node.score:.5f}) in {data[data >= node.score].count()}/{random_label_samples} cases.')
    print(f'Excpected p value from exponential distribution with u={dist.mean():.5f} and std={dist.std():.5f} p <= {(1-dist.cdf(node.score)):.5f}')
    
    if not save_intermediate_files:
        subprocess.run(f'rm {merged_fasta}* ', shell=True)
    
    return node, p, data, dist, tree, cropped_tree


# ---- OBSOLETE -----

# def tree_analysis(treefile, max_tree_leaves=1500,
#                             outlier_inv_gamma_low=0,
#                             outlier_inv_gamma_high=0.99,
#                             prok_min_size=2, prok_min_purity=0.5, euk_min_size=3, euk_min_purity=0.76):
#
#
#
#     # configure paths and flags
#     thread = current_process().pid
#     threadID_string = f'{thread} | :'
#
#     from ete3 import Tree
#     from core_functions.tree_functions import dirty_phyla_add, crop_leaves_to_size, get_outlier_nodes_by_invgamma, \
#         get_multiple_soft_LCAs
#     from paths_and_parameters import euk72_ep_tax, prok2111_as_tax
#
#     # load parsed taxonomy data
#     # standardize and reformat tax data
#     print(threadID_string + f' Parsing taxonomy from {prok2111_as_tax}')
#     prok_tax = pd.read_csv(prok2111_as_tax, sep='\t', names = ['acc', 'orgid', 'superkingdom', 'class'], index_col=0)
#
#     print(threadID_string + f' Parsing taxonomy from {euk72_ep_tax}')
#     euk_tax = pd.read_csv(euk72_ep_tax, sep='\t', names = ['acc', 'orgid', 'superkingdom', 'class'], index_col=0)
#
#     tax_merge = pd.concat([euk_tax, prok_tax])
#
#
#     # initialize tree
#     print(threadID_string + f' Reading Tree from {treefile}')
#     tree = Tree(treefile)
#
#     # name tree from file and header
#     # euk_header = load_pkl(euk72_header)
#     # tree_name = treefile.split('/')[1]
#     # tree_header = euk_header[euk_header.acc == tree_name].header.values
#
#     dirty_phyla_add(tree, tax_merge)
#
#     # merge leaf pairs until total amount of leaves is smaller than x
#     tree_size = len(tree.get_leaves())
#     if tree_size > max_tree_leaves:
#         print(threadID_string + f' Reducing tree size from {tree_size} to {max_tree_leaves}')
#         tree, crop_dict = crop_leaves_to_size(tree, max_tree_leaves)
#     else:
#         print(threadID_string + f' Tree of size {tree_size} has less than {max_tree_leaves} leaves, no cropping needed')
#
#     # calculate devaiting branch distances
#     print(threadID_string + f' Identifying outlier nodes by branch inverse gamma distribution')
#     outlier_nodes = get_outlier_nodes_by_invgamma(tree, p_low=outlier_inv_gamma_low, p_high=outlier_inv_gamma_high,
#                                                   only_leaves=False)
#
#     # cut_nodes = [node.detach() for node in outlier_nodes]
#
#     # calculate soft LCA nodes for prok and euk using partition entropy
#     print(threadID_string + f' Evaluating soft LCAs from {treefile}')
#
#     filter_taxa = set([leaf.tax_filter for leaf in tree.get_leaves()])
#     hard_LCA_dict = {}
#     soft_LCA_dict = {}
#
#     for tax in filter_taxa:
#         # hard_LCA_dict[tax] = get_paraphyletic_groups(tree, attribute='tax_filter', attr_value=tax)
#         soft_LCA_dict[tax] = get_multiple_soft_LCAs(tree, attribute='tax_filter', attr_value=tax,
#                                                     min_size=prok_min_size, min_purity=prok_min_purity)
#
#     valid_prok_LCAs = [key for key, value in soft_LCA_dict.items() if value != []]
#
#     soft_LCA_dict['Eukaryota'] = get_multiple_soft_LCAs(tree, attribute='tax_filter', attr_value='Eukaryota',
#                                                         min_size=euk_min_size, min_purity=euk_min_purity)
#
#     print(
#         threadID_string + f' Found valid soft LCAs for a total of {len(valid_prok_LCAs)} out of a possible {len(filter_taxa) - 1} taxons')
#     print(threadID_string + f' Found {len(soft_LCA_dict["Eukaryota"])} valid LCAs for Eukaryota')
#
#     # add LCA labels for LABEL visualisation
#     for node in tree.traverse():
#         node.add_feature('outlier', '')
#         node.add_feature('soft_LCA', '')
#         node.add_feature('soft_LCA_H', '')
#
#     # add the first valid soft_LCA as LCA for prok and all valid soft_LCAs for euk
#     for tax, nodes in soft_LCA_dict.items():
#         if nodes != []:
#             # add all euk nodes
#             if tax == 'Eukaryota':
#                 for node in nodes:
#                     node[0].soft_LCA = tax
#                     node[0].soft_LCA_H = node[3]
#
#                     # add first prok node
#             else:
#                 nodes[0][0].soft_LCA = tax
#                 nodes[0][0].soft_LCA_H = nodes[0][3]
#
#         # skip empty
#         else:
#             pass
#
#     #mark outliers
#     for node in outlier_nodes:
#         node.outlier = 'true'
#
#     # for more consistent visualisation
#     tree.ladderize()
#
#     # features [] is all features
#     tree.write(features=[], format_root_node=True,
#                outfile=treefile + '.annot')
#
#     return treefile + '.annot'
#
# # dirty tree printing
# def print_tree_to_pdf(treefile, view_in_notebook=False):
#
#     from ete3 import Tree, TreeStyle, NodeStyle, TextFace
#
#     tree = Tree(treefile, format=0)
#
#     # # add LCA nodes to list for NODE visualisation
#     # hard_LCA_nodes = [node for LCA_nodes in hard_LCA_dict.values() for node in LCA_nodes]
#     # # all first prok LCAs
#     # soft_LCA_nodes = [node[0][0] for node in soft_LCA_dict.values() if node != []]
#     # # add all euk LCAs
#     # soft_LCA_nodes = [node[0] for node in soft_LCA_dict['Eukaryota'] if node != []]
#
#     # define overall tree styling
#     ts = TreeStyle()
#     # ts.title.add_face(TextFace(tree_header, fsize=8), column=0)
#     ts.mode = 'r'
#     ts.show_leaf_name = False
#     ts.show_branch_length = False
#     ts.show_branch_support = True
#     # ts.optimal_scale_level = 'full'
#     ts.allow_face_overlap = True
#     ts.scale = 50
#
#     # set individual node styles
#     default_node_style = NodeStyle()
#     default_node_style['size'] = 0
#     default_node_style['fgcolor'] = 'Black'
#
#     default_leaf_style = NodeStyle()
#     default_leaf_style['size'] = 0
#     default_leaf_style['fgcolor'] = 'Black'
#
#     arc_node_style = NodeStyle()
#     arc_node_style['size'] = 2
#     arc_node_style['fgcolor'] = 'Cyan'
#
#     euk_node_style = NodeStyle()
#     euk_node_style['size'] = 2
#     euk_node_style['fgcolor'] = 'Green'
#
#     outlier_node_style = NodeStyle()
#     outlier_node_style['size'] = 3
#     outlier_node_style['fgcolor'] = 'Red'
#
#     hard_LCA_node_style = NodeStyle()
#     hard_LCA_node_style['size'] = 2
#     hard_LCA_node_style['fgcolor'] = 'Gray'
#
#     soft_LCA_node_style = NodeStyle()
#     soft_LCA_node_style['size'] = 4
#     soft_LCA_node_style['fgcolor'] = 'Black'
#
#     soft_euk_LCA_node_style = NodeStyle()
#     soft_euk_LCA_node_style['size'] = 4
#     soft_euk_LCA_node_style['fgcolor'] = 'Gray'
#
#     # set styling for all leaves and internal nodes
#     for node in tree.traverse():
#         node.set_style(default_node_style)
#
#         if node.is_leaf():
#             node.add_face(TextFace(node.name, fsize=4), column=1)
#             node.add_face(TextFace(' ' + node.tax_class, fsize=4), column=2)
#
#             if node.tax_superkingdom == 'Eukaryota':
#                 node.set_style(euk_node_style)
#
#             elif node.tax_superkingdom == 'Archaea':
#                 node.set_style(arc_node_style)
#
#         if node.soft_LCA != '':
#             node.add_face(TextFace(node.soft_LCA, fsize=4), column=0)
#             node.set_style(soft_LCA_node_style)
#
#         if node.outlier != '':
#             node.set_style(outlier_node_style)
#
#     # write to pdf
#     tree.render(treefile + '.pdf', tree_style=ts)
#
#     # for consistent view in notebook write to .png as well and view
#     if view_in_notebook:
#         from IPython.display import Image
#         if ts.mode == 'c':
#             tree.render(treefile + '.png', tree_style=ts, dpi=800, w=4000, h=4000)
#         else:
#             tree.render(treefile + '.png', tree_style=ts, dpi=800, w=1500, h=4000)
#
#         return Image(treefile + '.png')
#
#
# def realign_and_filter(query, root, threads=1, max_euk_leaf_size=200, max_prok_leaf_size=1300,
#                                  filter_entropy=0.5, muscle_reps_euk=25, muscle_reps_prok=5):
#
#
#     from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy
#     from core_functions.software_wrappers import muscle_ensamble
#
#     # configure paths and flags
#     thread = current_process().pid
#     threadID_string = f'{thread} | {query}:'
#
#     query_root = root + query + '/'
#     euk_fasta = query_root + query + '.fasta'
#     prok_fasta = query_root + query + '.members.fasta'
#
#     # for euk
#     euk_seqs = fasta_to_dict(file=euk_fasta)
#     euk_size = len(euk_seqs.keys())
#
#     # if there are too many eukaryotic sequences, crop to size
#     if euk_size > max_euk_leaf_size:
#         print(
#             threadID_string + f' There are more than {max_euk_leaf_size} sequences in {euk_fasta} ({euk_size}), will crop to size')
#         euk_fasta = fasta_reduce_size(euk_fasta, threads, max_euk_leaf_size, filter_entropy)
#
#     else:
#         print(
#             threadID_string + f' There are less than {max_euk_leaf_size} sequences in {euk_fasta} ({euk_size}), no cropping needed')
#
#         # align using muscle
#     print(threadID_string + f' Aligning with muscle5 as ensemble with {muscle_reps_euk} replicates')
#     euk_muscle = muscle_ensamble(euk_fasta, threads, muscle_reps_euk, super5=False)
#
#     # for prok
#     prok_seqs = fasta_to_dict(file=prok_fasta)
#     prok_size = len(prok_seqs.keys())
#
#     # if there are too many prokaryotic sequences, crop to size
#     if prok_size > max_prok_leaf_size:
#         print(
#             threadID_string + f' There are more than {max_prok_leaf_size} sequences in {prok_fasta}({prok_size}), will crop to size')
#         prok_fasta = fasta_reduce_size(prok_fasta, threads, max_prok_leaf_size, filter_entropy)
#
#     else:
#         print(
#             threadID_string + f' There are less than {max_prok_leaf_size} sequences in {prok_fasta}({prok_size}), no cropping needed')
#
#         # align using muscle
#     print(threadID_string + f' Aligning with muscle5 as ensemble with {muscle_reps_prok} replicates')
#     prok_muscle = muscle_ensamble(prok_fasta, threads, muscle_reps_prok, super5=False)
#
#     # filter both cropped alignments by columnwise bitscore
#     print(threadID_string + f' Filtering alignments to columnwise bitscore > {filter_entropy}')
#
#     # backup original alignments
#     subprocess.run(f'cp {euk_fasta}.muscle {euk_fasta}.muscle.b'.split())
#     subprocess.run(f'cp {prok_fasta}.muscle {prok_fasta}.muscle.b'.split())
#
#     aln = fasta_to_dict(file=f'{euk_fasta}.muscle')
#     aln_filter = filter_by_entropy(aln, filter_entropy)
#     dict_to_fasta(aln_filter, write_file=f'{euk_fasta}.muscle')
#
#     aln = fasta_to_dict(file=f'{prok_fasta}.muscle')
#     aln_filter = filter_by_entropy(aln, filter_entropy)
#     dict_to_fasta(aln_filter, write_file=f'{prok_fasta}.muscle')
#
#     return f'{euk_fasta}.muscle'
