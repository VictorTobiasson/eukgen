#storage of old and unwated functions


#several different weighting functions for attempting to evaluate cluster
#tree partitions using distortion and mututal information
#works fine but is currently not neccesary

# each leaf has a weight of its last node depth,
# sum weights of daughters postorder then normalize by total weight
# ---most consistent one so far---
def weight_tree_nodes_bottom_up(tree, add_residual=True):
    root = tree.get_tree_root()

    # add residual to correct for branch lengths of 0
    residual = 0
    if add_residual:
        residual = tree.get_farthest_node()[1] / len(tree)

    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.add_feature('weight', node.get_distance(root))
            # print(f'{node.name} has weight {node.weight}')
        else:
            weight = sum([child.weight + residual for child in node.children])
            node.add_feature('weight', weight)
            # print(f'Internal {node.name} has weight {node.weight}')

    total_weight = sum([child.weight for child in tree.children])

    for node in tree.traverse(strategy='postorder'):
        node.weight = node.weight / total_weight

    return total_weight


# leaves have weight of distance to root
# nodes have summed weight of children minus own distance to root
def weight_tree_nodes_bottom_up_corrected(tree, add_residual=True):
    root = tree.get_tree_root()

    # add residual to correct for branch lengths of 0
    residual = 0
    if add_residual:
        residual = tree.get_farthest_node()[1] / len(tree)

    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.add_feature('weight', node.get_distance(root) + residual)
            # print(f'{node.name} has weight {node.weight}')
        else:
            weight = node.get_distance(root) + sum([child.dist for child in node.children])
            node.add_feature('weight', weight + residual)
            # print(f'Internal {node.name} has weight {node.weight}')

    total_weight = sum([child.weight for child in tree.children])

    for node in tree.traverse(strategy='postorder'):
        node.weight = node.weight / total_weight

    return total_weight

# each node distributes weight among its decendants dependent on distance
def weight_tree_nodes(tree, add_residual=True):
    root = tree.get_tree_root()

    # add residual to correct for branch lengths of 0
    residual = 0
    if add_residual:
        residual = 1 / tree.get_farthest_node()[1] / len(tree)

    for node in tree.traverse(strategy='levelorder'):
        if node == root:
            # print('node is root weight is 1')
            tree.add_feature('weight', 1)

        # print(f'node has {len(node.children)} children')
        total_dist = sum([child.dist + residual for child in node.children])

        # print(total_dist)
        for child in node.children:
            weight = node.weight * child.dist / total_dist
            child.add_feature('weight', weight)


# leaves have weights of distance to root
# nodes have weights of sum of children
# then all nodes children get probabilities assigned as node.prob*(node.weight-child.weight)/node.weight
def weight_tree_nodes_up_down(tree, add_residual=True):
    root = tree.get_tree_root()
    total_weight = 0
    # add residual to correct for branch lengths of 0

    residual = 0
    if add_residual:
        residual = tree.get_farthest_node()[1] / len(tree)

    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.add_feature('weight', node.get_distance(root))
        else:
            weight = sum([child.weight + residual for child in node.children])
            node.add_feature('weight', weight)

    root.prob = 1

    for node in tree.traverse(strategy='levelorder'):
        # print(node.name, node.weight, node.prob)
        if len(node.children) == 1:
            # for linked internal nodes without branches
            child.add_feature('prob', node.prob)
        else:
            for child in node.children:
                prob = node.prob * (node.weight - child.weight) / node.weight
                child.add_feature('prob', prob)
                # print(child.name, child.prob)

    return total_weight


# leaves have weight equal to the distance to last node
# nodes have weights equal to sum of children
# normalised by total weight
def weight_tree_nodes_local(tree, add_residual=True):
    root = tree.get_tree_root()

    residual = 0
    if add_residual:
        residual = tree.get_farthest_node()[1] / len(tree)

    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.add_feature('weight', node.dist + residual)
        else:
            weight = sum([child.weight + residual for child in node.children])
            node.add_feature('weight', weight)

    for node in tree.traverse(strategy='levelorder'):
        sum_weight = sum([child.dist for child in node.children]) + node.dist
        for child in node.children:
            prob = child.dist / sum_weight
            child.add_feature('prob', prob)

        # leaves have weights equal to the root distance normalised of all leaves


# nodes have weights equal to the sum of their children
def weight_tree_nodes_top_down_prob(tree, add_residual=True):
    root = tree.get_tree_root()

    residual = 0
    if add_residual:
        residual = tree.get_farthest_node()[1] / len(tree)

    root.add_feature('weight', 1)

    for node in tree.traverse(strategy='levelorder'):
        sum_dists = sum([child.dist + residual for child in node.children])

        # for nested unbranched nodes in pathological cases
        if len(node.children) == 1:
            node.children[0].add_feature('weight', node.weight)

        else:
            for child in node.children:
                # farther nodes are lesss probable
                weight = node.weight * (sum_dists - child.dist + residual) / sum_dists
                # farther nodes are more probable
                # weight = node.weight*(child.dist+residual)/sum_dists
                child.add_feature('weight', weight)


# using precomputed weights calculate the MI and distortion
#meant for automated cluster evalutaion on trees
def calculate_mututal_I_and_distortion(tree, nodes, uniform=False):
    # initialise
    leaves = tree.get_leaves()
    depth = tree.get_farthest_leaf()[1]
    Icx = 0
    dis = 0
    tiny = np.finfo(np.float64).eps

    Px = np.array([leaf.weight for leaf in leaves])

    for C in nodes:
        Pc = C.weight
        # mututal information is done over all x
        Plca2 = np.array([C.get_common_ancestor(x).weight ** 2 for x in leaves])
        Icx += (((Pc ** 2) / Plca2) * Px * np.log2(Pc / Plca2)).sum()

        # distortion calculations are done over all x
        # difference in probabilities
        d = np.abs(Px - Pc)
        dis += (((Pc ** 2) / Plca2) * Px * d).sum()

    return dis, Icx



# remove isolated leaves by checking all monophyletic nodes for singletons
#currently unused way of deleting singletons leaves
#superceded by "get_soft_LCA_by_relative_entropy"
def trim_singleton_leaves(tree, attribute='tax_superkingdom', attr_value='Eukaryota', min_size=1, detach=True):
    LCA_groups = get_paraphyletic_groups(tree, attribute=attribute, attr_value=attr_value)

    pruned_leaves = []

    for node in LCA_groups:
        # remove all LCA leaves which are not in the majority among its neighbors
        if len(node.get_leaves()) <= min_size:
            neigbour_attr = [getattr(leaf, attribute) for leaf in node.up.get_leaves()]

            if neigbour_attr.count(attr_value) / len(neigbour_attr) < 0.5 and detach:
                pruned_leaves.append(node.detach())

            else:
                pruned_leaves.append(node)

    return pruned_leaves



#from ipynb, meant to try to evaluate taxonomic distribution across seqrch hits without constructing trees.
#Does is not sensitive enough, hence the full tree construction or microscosm pipeline
# given a series of query accessions, retreive all proteins from prokaryotic hits and return their phylogenetic distribution count
# designed for multiprocess Pool and checkpointing to .pkl files
# not very portable

def get_hit_tax_dist(queries):
    stime = time.time()
    querynr = len(queries)

    thread = multiprocessing.current_process().pid
    print(f'{thread}: started\n')

    print(f'{thread}: processing {querynr} queries\n')

    # load data from files
    print(f'{thread}: loading search\n')
    searchDF = load_pkl(root + 'analysis/core_data/merged_filtered_cov20_self-match_tsv_edited_no_aln.pkl')

    print(f'{thread}: loading clust\n')
    prok_clust = load_pkl(root + 'analysis/core_data/prok2111_filtered-prof-search-clust.pkl')['members']
    reindex(prok_clust, 'cluster_acc')

    print(f'{thread}: loading tax\n')
    prok_tax = load_pkl('analysis/core_data/prok2111_protein_taxonomy_trimmed.pkl')

    # process only the given slice
    searchDF = searchDF[(searchDF.Query.isin(queries)) &
                        (searchDF.Query != searchDF.Target) &
                        (searchDF.Pairwise_cov > 0.5)]

    # set index for faster iterating over queries
    searchDF.sort_values(by='Query', inplace=True)
    searchDF.set_index(keys=['Query'], drop=True, inplace=True)

    # iterate and pool taxa distributions
    taxa = {}
    n = 0
    printn = 50
    checkn = 1000

    for query in queries:

        if n % printn == 0:
            print(f'{thread}: calculating {query} \t{n}|{querynr} \tT+{round(time.time() - stime)} seconds\n')

        # find all target profile hits
        profiles = pd.Series(searchDF.loc[query, 'Target'])

        # find all proteins in target profile hits
        proteins = pd.Series(prok_clust.loc[profiles, 'acc'])

        # find all taxonomic information from proteins in taget profile hits
        query_taxa = prok_tax.loc[proteins, 'class'].value_counts()

        # add to dict
        taxa[query] = query_taxa

        # itermediate save
        if n != 0 and n % checkn == 0:
            print(f'{thread}: saved checkpoint {n / checkn}')
            dump_pkl(taxa, f'analysis/core_data/tax/{thread}_checkpoint_{int(n / checkn)}_tax.pkl')
            taxa = {}

        n += 1

    dump_pkl(taxa, f'analysis/core_data/tax/{thread}_tax.pkl')
    return taxa


# launch parallel execution
queries = searchDF[(searchDF.Query != searchDF.Target)
                   & (searchDF.Pairwise_cov > 0.5)].Query.unique()

splits = np.array_split(queries, 16)
with multiprocessing.Pool(processes=16) as pool:
    pool.map(get_hit_tax_dist, splits)

# load data from savepoints into one dictionary
tax_data = {}
for file in os.listdir('analysis/core_data/tax/'):
    print('analysis/core_data/tax/' + file)
    data = load_pkl('analysis/core_data/tax/' + file)
    tax_data = tax_data.copy()
    tax_data.update(data)

# merge dict series into one dataframe
query_tax = pd.DataFrame()
for query, data in tax_data.items():
    print(query)

    data.name = query

    temp_tax = pd.DataFrame(data).transpose()
    query_tax = pd.concat([query_tax, temp_tax])

# save processed dataframe
dump_pkl(query_tax, 'analysis/core_data/hit_distribution_cov50.pkl')

# calculate query statistics tables
# load raw hit counts
query_tax = load_pkl('analysis/core_data/hit_distribution_cov50.pkl')

# normalize to relative total hits
query_tax_rel = query_tax.div(query_tax.sum(axis=1), axis=0)
# calculate percentile ranks for relative observations skipping 0 observations
query_tax_rel_percentile = query_tax_rel.apply(lambda df: df[df != 0].rank(method='max', pct=True), axis=0).fillna(0)

# multiply by individual hits to get weights
query_tax_weight = pd.DataFrame(query_tax.values * query_tax_rel.values, columns=query_tax.columns,
                                index=query_tax.index)
# calculate percentile rank of observation of weight
query_tax_weight_percentile = query_tax_weight.apply(lambda df: df.rank(method='max', pct=True), axis=0)

# write proteins from query clusters which have top 10% relative hits in
# respective taxon to file

prominent_taxon_proteins = {}
for taxon in query_tax.columns:
    prominent_queries = query_tax_rel_percentile[(query_tax_rel_percentile[taxon].between(0.85, 0.90))].index
    prominent_proteins = euk_clust.loc[prominent_queries, 'acc'].values

    homo_tax = euk_tax[euk_tax['class'] == 'Mammalia']
    prominent_homo_proteins = homo_tax[homo_tax.index.isin(prominent_proteins)].index
    outdata = euk_header[euk_header.acc.isin(prominent_homo_proteins)]
    print(taxon, outdata.shape[0])
    outdata.reset_index(inplace=True)

    outdata['header'].to_csv(f'analysis/core_data/significant/{taxon.replace(" ", "_")}.tsv', sep='\t', index=None,
                             header=None)

