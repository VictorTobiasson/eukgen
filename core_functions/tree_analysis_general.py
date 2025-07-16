# temporary cope of tree analysis function for tree annotation using a ARCHAEAL OUTGROUP
# SUBSTITUTING EUKARYA FOR OUTGROUP
# assigns taxonomy to leaves, identifies and trims outliers, assigns LCA nodes and performs distance calculations
# outputs a treeDF with tabulated data for each eukaryotic clade and its corresponding prokaryotic sister-clades
def tmp_tree_analysis_general(tree_file, leaf_mapping, tree_name, outgroup,
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
    from core_functions.microcosm_functions import retrieve_LCA_scope

    # load tree
    print(f'Reading Tree from {tree_file}')
    tree = Tree(tree_file)

    # load leaf mapping from tree cropping
    print(f'Loading leaf data from {leaf_mapping}')
    leafDF = pd.read_csv(leaf_mapping, sep='\t', names=['acc', 'leaf', 'superkingdom', 'full_class', 'rank'],
                         index_col=0)

    # merge all euk classification in to 'Eukaryota' save old class
    leafDF['class'] = leafDF.full_class
    # leafDF.loc[leafDF[leafDF.superkingdom == 'Eukaryota'].index, 'class'] = 'Eukaryota'

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
    taxa_list = leafDF[(leafDF.superkingdom.isin(['Bacteria', 'Archaea'])) 
                       & (leafDF.leaf != 'DELETED')
                       & (leafDF['class'] != outgroup)
                      ]['class'].unique()

    lca_dict = {}
    for taxa in sorted(taxa_list):
        lca_dict[taxa] = get_multiple_soft_LCA_by_relative_purity(tree, taxa,
                                                                  n_best=9999,
                                                                  min_size=prok_clade_size,
                                                                  min_purity=prok_clade_purity,
                                                                  max_depth=10)

    # find all euk sof LCAs meeting criteria
    lca_dict[outgroup] = get_multiple_soft_LCA_by_relative_purity(tree, outgroup,
                                                                     n_best=9999,
                                                                     min_size=euk_clade_size,
                                                                     min_purity=euk_clade_purity,
                                                                     max_depth=10)

    # display LCA results
    accepted_taxa = [taxa for taxa, node in lca_dict.items() if node != [] and taxa != outgroup]
    rejected_taxa = [taxa for taxa, node in lca_dict.items() if node == [] and taxa != outgroup]

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

    print(f'\nFound {len(lca_dict[outgroup])} LCA nodes for outgroup')
    for node_data in lca_dict[outgroup]:
        print(f'    {outgroup:<25} size: {node_data[1]}\t weight: {round(node_data[3], 2)}')
        pass
    print()

    # if more than three good euk nodes are identified warn and restruct analysis
    if len(lca_dict[outgroup]) > 3:
        print(
            f'WARNING! High paraphyly in {outgroup}, considering the three largest clades of {len(lca_dict[outgroup])} total!')
        lca_dict[outgroup] = lca_dict[outgroup][:3]

    # revert to original eukaryotic classes on tree
    leafDF.columns = ['leaf', 'superkingdom', 'class', 'rank', 'filter_class']
    tree = map_leafDF(tree, leafDF)

    # calculate distances and format DataFrame
    euk_dist_dict = {}

    # separate euk and prok LCA nodes
    all_LCA_data = [node for taxa, nodes in lca_dict.items() for node in nodes]
    all_LCA_nodes = [node[0] for node in all_LCA_data]
    euk_LCA_nodes = [node for node in all_LCA_data if node[0].lca_taxa == outgroup]
    prok_LCA_nodes = [node for node in all_LCA_data if node[0].lca_taxa != outgroup]

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
        f'Consider only closest clades for downstream analysis, clade limit per {outgroup} clade is {consider_closest_n_prok_LCAs}')

    for euk_clade_rep, data in tree_data.groupby('euk_clade_rep'):
        print(f'Including the following LCA list for {outgroup} LCA: {euk_clade_rep}')
        data.iloc[0].euk_node_ref.add_feature('LCA', outgroup)

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


# generate a set of constraint trees for each euk_LCA * prok_LCA pair for IQtree2 -g constraint testing
# ((euk_LCA, prok_sister),(remaining_prok));
# return a dataframe of jobs to be run
# clade size dictates the maximum number of sequences forming the contrained clades
# CURRENLTY ONLY ALLOWING FOR VARIABLE EUK CLADE SIZE ONLY USING ONE PROK SEQUENCE
def tmp_format_constraint_analysis_general(root, basename, tree_data, outgroup, clade_size=10):

    import subprocess
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
            if node.LCA == outgroup:
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

    print(f'Formatted {constraint_job_data.shape[0]} guide trees for {len(euk_LCA_nodes)} {outgroup} LCAs')

    return constraint_job_data.set_index('euk_clade_rep')
