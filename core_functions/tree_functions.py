
import pandas as pd
import numpy as np
from scipy import stats


# TREE OPERATIONS

# fit inverse gamma distribution to all internal node distances and
# exclude those beyond threshold probability
def get_outlier_nodes_by_invgamma(tree, p_low=0, p_high=0.99, only_leaves=False, deletion_cutoff=0.15):
    from scipy import stats

    tree_size = len(tree)

    if only_leaves:
        node_dists = [(node, node.dist) for node in tree.get_leaves()]

    else:
        node_dists = [(node, node.dist) for node in tree.traverse()]

    dist_series = pd.Series([i[1] for i in node_dists])

    fit_alpha, fit_loc, fit_beta = stats.invgamma.fit(dist_series.values)

    cutoff_high = stats.invgamma.ppf(p_high, a=fit_alpha, loc=fit_loc, scale=fit_beta)
    cutoff_low = stats.invgamma.ppf(p_low, a=fit_alpha, loc=fit_loc, scale=fit_beta)

    outlier_nodes = [node[0] for node in node_dists if node[1] < cutoff_low or node[1] > cutoff_high]
    print(f'Identified {len(outlier_nodes)} outlier nodes outside interval {cutoff_low} > d > {cutoff_high}')

    bad_leaves = set()
    for node in outlier_nodes:
        bad_leaves = bad_leaves.union({*[leaf for leaf in node.get_leaves()]})

    if len(bad_leaves) > tree_size * deletion_cutoff:
        print(
            f'WARNING: Outlier nodes map to more than {deletion_cutoff * 100}% ({len(bad_leaves) / tree_size * 100}%) of the tree leaves, refusing to crop')
        outlier_nodes = []

    return outlier_nodes, dist_series, (fit_alpha, fit_loc, fit_beta), (cutoff_high, cutoff_low)


# fit lognorm distribution to all internal node distances, add small branch length to avoid log(0)
# exclude those beyond threshold probability
def get_outlier_nodes_by_lognorm(tree, p_low=0, p_high=0.99, only_leaves=False, deletion_cutoff=0.15, delete=True):
    from scipy import stats

    tree_size = len(tree)

    if only_leaves:
        node_dists = [(node, node.dist + 0.00000001) for node in tree.get_leaves()]

    else:
        node_dists = [(node, node.dist + 0.00000001) for node in tree.traverse()]

    dist_series = pd.Series([i[1] for i in node_dists])

    fit_alpha, fit_loc, fit_beta = stats.lognorm.fit(dist_series.values)

    cutoff_high = stats.lognorm.ppf(p_high, s=fit_alpha, loc=fit_loc, scale=fit_beta)
    cutoff_low = stats.lognorm.ppf(p_low, s=fit_alpha, loc=fit_loc, scale=fit_beta)

    outlier_nodes = [node[0] for node in node_dists if node[1] < cutoff_low or node[1] > cutoff_high]
    print(f'Identified {len(outlier_nodes)} outlier nodes outside interval {cutoff_low} > d > {cutoff_high}')

    bad_leaves = set()
    for node in outlier_nodes:
        bad_leaves = bad_leaves.union({*[leaf for leaf in node.get_leaves()]})

    if len(bad_leaves) > tree_size * deletion_cutoff:
        print(
            f'WARNING: Outlier nodes map to more than {deletion_cutoff * 100}% ({len(bad_leaves) / tree_size * 100}%) of the tree leaves, refusing to crop')
        outlier_nodes = []

    return outlier_nodes, dist_series, (fit_alpha, fit_loc, fit_beta), (cutoff_high, cutoff_low)

# reroot tree inplace so that each side of the root has equal stem length mass
def weighted_midpoint_root(tree):
    # reroot at midpoint as starting reference
    root_point = tree.get_midpoint_outgroup()
    tree.set_outgroup(root_point)

    # precalculate decendant sum of dists
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            # print(f'node is leaf adding w {node.dist}')
            node.add_feature('child_weight', node.dist)

        else:
            child_weight = sum([child.child_weight for child in node.children])
            # print(f'node is not leaf w {child_weight}')
            node.add_feature('child_weight', child_weight)

    total_weight = tree.child_weight
    best_balance = total_weight
    best_node = None

    tree.add_feature('balance', total_weight)

    # start by adding tree root children
    nodes_to_consider = [child for child in tree.children]

    # while there are nodes to consider
    while nodes_to_consider:

        node = nodes_to_consider.pop()
        my_weight = node.child_weight

        # parent_weight is total_weight - my_weight, balance is parent_weight - my_weight
        balance = abs(total_weight - 2 * my_weight)
        node.add_feature('balance', balance)

        # if my balance is better than my parents consider my children
        if balance <= node.up.balance:
            nodes_to_consider.extend(node.children)

        # if current node is best so far update best node
        if balance <= best_balance:
            best_balance = balance
            best_node = node

    # set outgroup and adjust branch lengths
    tree.set_outgroup(best_node)

    # calculation of centerpoint
    # dist1+dist2 = length
    # weight1*dist1 = weight2*dist2
    # dist1 = weight2*Length/(weight2+weight1)
    # dist2 = weight1*Length/(weight2+weight1)

    sum_length = sum([node.dist for node in tree.children])
    sum_balance = sum([node.child_weight for node in tree.children])
    for node in tree.children:
        node.dist = sum_length - (node.child_weight * sum_length / sum_balance)

    return tree

# consider all "cherry nodes" with exactly two leaf children
# iteratively crop most similar cherry clades until size
# monophyletic=True requires tree feature ["tax"] for every leaf
# if monophyletic only cherries with same tax are considered
def crop_leaves_to_size(tree, max_size, save_cropped_references=True, monophyletic=False):
    # helper function to assess whether node is a cherry
    def node_is_cherry_root(cherry_root, monophyletic=False):
        cherry_children = cherry_root.children

        # if the node has exactly two leaf children its a "cherry"
        if [node.is_leaf() for node in cherry_children] == [True, True]:

            # if running as monophyletic and children do not share taxa skip
            if monophyletic and cherry_children[0].taxa != cherry_children[1].taxa:
                return

            # format the cherry and return
            cherry = (cherry_root,
                      cherry_children[0].get_distance(cherry_root),
                      cherry_children[1].get_distance(cherry_root))

            return cherry

        # retrun None if not cherry root
        else:
            return

            # calculate number of leaves to remove to reach max_size

    leaves = tree.get_leaves()
    leaves_to_remove = len(leaves) - max_size

    crop_dict = {leaf.name: [leaf.name] for leaf in leaves}

    # if tree is small enough
    if leaves_to_remove < 1:
        return tree, crop_dict

    # find and mark all initial cherries (nodes with exactly two leaf children)
    cherries = []
    while leaves != []:

        # consider nodes in reverse order for speed
        cherry_root = leaves.pop().up

        # if the node has exactly two leaf children its a "cherry"
        cherry = node_is_cherry_root(cherry_root, monophyletic)

        if cherry != None:
            cherries.append(cherry)

            # second cherry leaf should not be considered, pop it as well.
            leaves.pop()

    # sort list of initial cherries
    cherries = sorted(cherries, key=lambda x: x[1] + x[2], reverse=True)

    # repeat until cropped
    while leaves_to_remove > 0:

        if cherries == []:
            return tree, crop_dict

        cherry = cherries.pop()

        # delete farthest leaf
        if cherry[1] > cherry[2]:
            crop_leaf = cherry[0].get_children()[0]
            save_leaf = cherry[0].get_children()[1]

        else:
            crop_leaf = cherry[0].get_children()[1]
            save_leaf = cherry[0].get_children()[0]

        # the retained leaf remembers the cropped leaf name and all its deleted leaves
        if save_cropped_references:
            crop_dict[save_leaf.name] = crop_dict[save_leaf.name] + crop_dict[crop_leaf.name]

        # remove cropped leaf and its dict entry
        crop_dict.pop(crop_leaf.name, None)
        crop_leaf.delete()

        # deletion of the leaf might have created a new cherry to be conisdered
        cherry_root = save_leaf.up
        cherry = node_is_cherry_root(cherry_root, monophyletic)

        # if a new cherry is created insert the cherry into the sorted list
        if cherry != None:
            cherry_distance = cherry[1] + cherry[2]


            # check each element if the cherry distance is larger than for the new cherry
            # start from smaller values as this is likley faster
            # this runs in O(N) time but could be improved to N(log(N)) with binary insertion
            # python bisect module does this but 3.9 implementation does not take key= so its annoying
            insert_i = 0
            for i in reversed(range(len(cherries))):
                i_distance = cherries[i][1] + cherries[i][2]
                if i_distance > cherry_distance:
                    insert_i = i
                    break

            # insert
            cherries = cherries[:insert_i + 1] + [cherry] + cherries[insert_i + 1:]

        leaves_to_remove -= 1

    return tree, crop_dict


# take dict of leaf_acc:[collapsed_leaf, collapsed_leaf, collapsed_leaf ...]
# transform into longfor annotated DF for tax mapping
# taxDF is acc: [orgid, superkingdom, class]
def format_leafDF(crop_dict, taxDF):
    # convert crop_dict to longform DF
    leafDF = pd.DataFrame([[leaf, acc] for leaf, accs in crop_dict.items() for acc in accs if accs], columns=['leaf', 'acc'])

    leafDF = leafDF.set_index('acc')

    # construct leaf to tax mapping
    leafDF[['superkingdom', 'class']] = taxDF.loc[leafDF.index][['superkingdom', 'class']]

    return leafDF


# update existing leafDF with new dict of cropped leaves
def update_leafDF(crop_dict, leafDF):
    inverse_dict = {acc: leaf for leaf, accs in crop_dict.items() for acc in accs}
    leafDF['leaf'] = [inverse_dict[acc] for acc in leafDF['leaf']]

    return leafDF


# map data from leafDF to tree leaves
def map_leafDF(tree, leafDF):
    tree_leaf_dict = {leaf.name: leaf for leaf in tree.get_leaves()}

    leafDF = leafDF.reset_index().set_index('leaf')

    for leaf_acc in tree_leaf_dict.keys():
        # reformat value_counts as explicit pd.series to avoid string return for single entries
        tax_counts = pd.Series(leafDF.loc[leaf_acc]['class']).value_counts()

        tree_leaf_dict[leaf_acc].add_feature('taxa', list(tax_counts.index.values))
        tree_leaf_dict[leaf_acc].add_feature('counts', list(tax_counts.values))

    return tree


# replace an LCA node with its closes leaf, tracking mergers in crop_dict
def collapse_LCA(node, crop_dict):
    # find closest LCA leaf node to save
    # could be replaced with other criteria
    save_leaf = node.get_closest_leaf()[0]

    # list all LCA leaves which are to be removed
    leaf_names = [leaf.name for leaf in node.get_leaves() if leaf.name != save_leaf.name]

    # extend saved crop_dict entry with entries from all removed leaves
    deleted_leaves = []
    for name in leaf_names:
        crop_dict[save_leaf.name].extend(crop_dict.pop(name))
        # crop_dict[save_leaf.name].extend(crop_dict.pop(name, None))

    # delete bottom up to pererve various topologies
    for leaf in node.get_leaves():
        if leaf != save_leaf:
            leaf.delete()

    return crop_dict


# requires annotated tree with leaf.taxa with list of taxa mapping for collapsed leaf,
# and leaf.counts with list of frequencies of those taxa
def get_multiple_soft_LCA_by_relative_purity(tree, tax, n_best, min_size, min_purity, max_depth=-1,
                                             verbose_leaf_stats=False):
    # inital values
    all_leaves = tree.get_leaves()
    nodes_found = 0
    soft_LCA_nodes = []

    # mark all leaves as unchecked
    for leaf in all_leaves:
        leaf.add_feature('to_check', True)

    # consider all points from leaf with atleast one sequence in tax to root if not already checked
    tax_leaves = [leaf for leaf in all_leaves if tax in leaf.taxa and leaf.to_check]

    while nodes_found < n_best:

        # of no more leaves found
        if len(tax_leaves) == 0:
            return soft_LCA_nodes

        # get set of all nodes on a path from any tax_leaf to root within depth limit
        check_nodes = {node for leaf in tax_leaves for node in leaf.get_ancestors()[:max_depth]}

        # add all leaves as they could contain collapsed clades
        check_nodes = check_nodes.union({leaf for leaf in tax_leaves if sum(leaf.counts) > 1})

        # calculate total unchecked labels on tree
        size_tree = sum([sum(leaf.counts) for leaf in all_leaves if leaf.to_check])

        # total unchecked labels of tax on tree
        size_tax = sum(
            [leaf.counts[leaf.taxa.index(tax)] if tax in leaf.taxa and leaf.to_check else 0 for leaf in tax_leaves])

        node_weight_data = []

        # calculate weight as (fraction of tax in clade)*(fraction of tax on tree) or purity*scope
        # single loop for speed
        for node in check_nodes:

            node_leaves = node.get_leaves()
            # total unchecked size of clade
            size_clade = 0

            # total number of unchecked tax in clade
            size_tax_in_clade = 0

            for leaf in node_leaves:
                size_clade += sum(leaf.counts)
                if tax in leaf.taxa:
                    size_tax_in_clade += leaf.counts[leaf.taxa.index(tax)]

            node_weight = (size_tax_in_clade / size_clade) * (size_tax_in_clade / size_tax)

            # save weight values if they are pure and large enough and weight is reasonable (rare cases of nested clades?)
            if (size_tax_in_clade / size_clade >= min_purity) and (size_tax_in_clade >= min_size) and (node_weight <= 1):
                node_weight_data.append([node, size_clade, size_tax_in_clade, node_weight])

        # best node is node with highest weight
        node_weight_data = sorted(node_weight_data, key=lambda x: x[-1], reverse=True)

        # if any nodes passed filters
        if node_weight_data != []:

            optimal_node = node_weight_data[0]

            # mark all leaves from optimal node as checked, not to considered in further rounds
            for leaf in optimal_node[0].get_leaves():
                leaf.to_check = False

            # update nodes to consider
            tax_leaves = [leaf for leaf in tax_leaves if tax in leaf.taxa and leaf.to_check]

            # count node and add to list of found nodes
            nodes_found += 1
            soft_LCA_nodes.append(optimal_node)

            if verbose_leaf_stats:
                print(f'Found best node for {tax}')
                print('Tree size:                              ', size_tree)
                print('Total tax:                              ', size_tax)
                print('Current clade size:                     ', optimal_node[1])
                print('Total tax in clade:                     ', optimal_node[2])
                print('Fraction of total tax covered by clade: ', optimal_node[2] / size_tax)
                print('Purity of clade:                        ', optimal_node[2] / optimal_node[1])
                print('Final weight:                           ', optimal_node[3])
                # print(optimal_node[0].get_ascii(attributes=['taxa','counts']))
                print()

        else:
            return soft_LCA_nodes

    return soft_LCA_nodes


# phylogenetically "aware" tree cropping procedure
# assignes taxa from taxDF to tree
# collpses monophyletic groups and then finds LCAs of given purity from remaining leaves
# returns a tree with the highest ranked clades by size
def crop_leaves_to_size_considering_taxa(tree, taxDF, max_size, min_clade_size=2, min_clade_purity=0.9,
                                         LCA_search_depth=3):
    # create initial naive leaf mapping
    crop_dict = {name: [name] for name in tree.get_leaf_names()}

    print(
        f'Reducing tree with {len(crop_dict.keys())} leaves to {max_size} maintaining a clade purity of {min_clade_purity}')

    leafDF = format_leafDF(crop_dict, taxDF)
    tree = map_leafDF(tree, leafDF)
    taxa_list = leafDF['class'].unique()

    # collapse all monophyletic clades
    tree, crop_dict = crop_leaves_to_size(tree, max_size, save_cropped_references=True, monophyletic=True)

    print(f'Enforcing monophyly across {len(taxa_list)} taxa reduced tree to {len(crop_dict.keys())} leaves')

    # update leafDF to remap leaves
    leafDF = update_leafDF(crop_dict, leafDF)
    tree = map_leafDF(tree, leafDF)

    taxa_list = leafDF['class'].unique()

    # collapse all identified taxa to purity limit
    n_best = 99999
    for taxa in taxa_list:
        lca_nodes = get_multiple_soft_LCA_by_relative_purity(tree, taxa, n_best, min_clade_size, min_clade_purity,
                                                             LCA_search_depth)

        # reduce LCA to representatives while maintaining leaf tax count
        for node_data in lca_nodes:

            lca_node = node_data[0]

            # skip LCAs which map to leaves as these cannot be reduced further
            # in edge cases LCAs can become roots following deleted nodes
            if lca_node.is_root() or lca_node.is_leaf():
                continue

            else:
                # delete leaves and track collapsed leaves
                crop_dict = collapse_LCA(lca_node, crop_dict)

    print(f'After collapsing soft LCAs of high purity tree now has {len(crop_dict.keys())} leaves')

    # update leafDF and remap counts to tree
    leafDF = update_leafDF(crop_dict, leafDF)
    tree = map_leafDF(tree, leafDF)

    # calculate ranks of all taxa sizes
    all_ranks = pd.DataFrame()
    for taxa in leafDF['class'].unique():
        tax_series = pd.DataFrame(leafDF[leafDF['class'] == taxa].leaf.value_counts())
        tax_series['rank'] = tax_series.rank(method='dense', pct=True)
        tax_series['class'] = [taxa for _ in tax_series.index]
        all_ranks = pd.concat([all_ranks, tax_series])

    # expand to full data series sorted by leaf acc
    rank_list = []
    all_ranks = all_ranks.reset_index().sort_values(by=['leaf', 'class'])
    for i, row in all_ranks.iterrows():
        rank_list.extend([row['rank']] * row['count'])

    # add to leafDF sorted by leaf acc
    leafDF = leafDF.sort_values(by=['leaf', 'class'])
    leafDF['rank'] = rank_list

    # good leaves are those which are among the top 500 scoring for the best taxa
    good_leaves = all_ranks.sort_values(by='rank', ascending=False).drop_duplicates(subset='leaf', keep='first')[0:max_size-1].leaf.values

    good_seqs = leafDF[leafDF.leaf.isin(good_leaves)].shape[0]

    print(f'After keeping the {max_size} leaves with highest taxa-specific size rank tree clades map to {good_seqs} original sequences {round(good_seqs / leafDF.shape[0] * 100, 1)}% of original tree')

    leafDF.leaf = [leaf if leaf in good_leaves else 'DELETED' for leaf in leafDF.leaf]

    for leaf in tree.get_leaves():
        if leaf.name not in good_leaves:
            leaf.delete()

    return tree, leafDF

# --------- OBSOLETE ----------

# write pdf tree under /data/tobiassonva/data/eukgen/tmp_trees
# overwrites previous treee in directory unless specified name
# the file path is always relative to the notbook starting directory due to the "Tornado" viewer restrictions
def view_tree_pdf(treefile):

    from IPython.display import IFrame
    from paths_and_parameters import path_notebooks
    import subprocess

    tmp_tree_file = path_notebooks + 'tmp_trees/tmp_tree.pdf'

    subprocess.run(f'ln -sf {treefile} {tmp_tree_file}'.split())

    IFrame(tmp_tree_file, width=950, height=950 * 1.5)

# return LCA for list of leaf names
def get_LCA(tree, leaf_names):
    leaves = [tree.get_leaves_by_name(name)[0] for name in leaf_names]
    LCA = leaves[0].get_common_ancestor(leaves)
    return LCA

#extract all against all distance matrix for ete3 trees
def calculate_pairwise_distances(tree):
    leaves = tree.get_leaves()
    pairwise_mat = np.zeros((len(leaves),len(leaves)))
    for i, m in enumerate(leaves):
        for j, k in enumerate(leaves[i+1:]):
            pairwise_mat[i,j] = tree.get_distance(m, k)
    return pairwise_mat


# return the set of all leaf pairs
def get_sister_leaf_sets(tree):
    leaves = set(tree.get_leaves())
    sister_sets = []
    while leaves:
        leaf = leaves.pop()
        sister = leaf.get_sisters()[0]
        if sister in leaves:
            sister_sets.append(set((leaf, sister)))

    return sister_sets


# give each node in tree a feature numerical id
def add_node_ids(tree, strategy='postorder'):
    for i, node in enumerate(tree.traverse(strategy=strategy)):
        node.add_feature('post_i', i)
    return tree



# count all combintions of nodes such that no node is a decendant of any other
# infeasible for more than 30 nodes as combinations scale fast!
def enumerate_all_clades(tree):
    stack = {}
    clades = []

    # stack a copy of the input tree as first tree
    # use ID to prevent duplicates
    subtree = tree.copy('deepcopy')
    stack[subtree.get_topology_id()] = subtree

    while stack:
        print(f'Stack contains {len(stack)}')

        # take the last element out of the stack by id
        subtree_id = list(stack.keys()).pop()
        subtree = stack.pop(subtree_id)

        # save the leaf configuration
        clades.append(tuple(leaf.post_i for leaf in subtree.get_leaves()))

        sister_sets = get_sister_leaf_sets(subtree)
        for pair in sister_sets:

            ids = [i.post_i for i in pair]

            recursetree = subtree.copy()

            # collapse sister nodes
            for node in recursetree.traverse():
                if node.post_i in ids:
                    node.detach()

            # add cropped tree to stack if the root has children
            if recursetree.children:
                stack[recursetree.get_topology_id()] = recursetree

    return clades


# iteratively merge closest leaf pair until less than N leaves
def reduce_leaves_to_size(tree, max_size):
    current_size = len(tree.get_leaves())

    if max_size >= current_size:
        print(f'Tree of length {current_size} smaller than {max_size}')
        return tree

    leaf_partners = {leaf.name: get_closest_leaf(leaf) for leaf in tree.get_leaves()}
    leaf_partner_dist = {key: value[1] for key, value in leaf_partners.items()}
    leaf_partners = {key: value[0].name for key, value in leaf_partners.items()}

    # print(leaf_partner_dist)
    # print(leaf_partners)

    # serial implementetaion, not ideal as it produced uneven pruning if terminal brach length are very even.
    while max_size < current_size:

        min_distance = min(leaf_partner_dist.values())
        min_leaf_A = list(leaf_partner_dist.keys())[list(leaf_partner_dist.values()).index(min_distance)]

        min_leaf_B = leaf_partners[min_leaf_A]

        # print(f'{current_size} checking {min_leaf_A}, deleting closest partner is {min_leaf_B} with distance {min_distance}')

        # delete closest leaf
        try:
            tree.get_leaves_by_name(min_leaf_B)[0].delete()

        # if current leaf A maps to a deleted leaf update closest leaf and delete
        except IndexError:
            # update the min_leaf with new closest pair
            new_leaf = tree.get_leaves_by_name(min_leaf_A)[0]
            closest_new_leaf = get_closest_leaf(new_leaf)

            leaf_partner_dist[new_leaf.name] = closest_new_leaf[1]
            leaf_partners[new_leaf.name] = closest_new_leaf[0].name
            min_leaf_B = leaf_partners[min_leaf_A]

            tree.get_leaves_by_name(min_leaf_B)[0].delete()

        # delete the removed partner from dictionaries
        leaf_partner_dist.pop(min_leaf_B, None)
        leaf_partners.pop(min_leaf_B, None)

        # update the min_leaf with new closest pair
        new_leaf = tree.get_leaves_by_name(min_leaf_A)[0]
        closest_new_leaf = get_closest_leaf(new_leaf)

        leaf_partner_dist[new_leaf.name] = closest_new_leaf[1]
        leaf_partners[new_leaf.name] = closest_new_leaf[0].name

        current_size -= 1

        # print(leaf_partner_dist)
        # print(leaf_partners)

    return tree


# quick function for adding phylogenetic annotation to tree labels
def dirty_phyla_add(tree, tax_mapping):
    for leaf in tree.get_leaves():
        leaf.add_feature('tax_superkingdom', 'NONE')
        leaf.add_feature('tax_filter', 'NONE')
        try:
            acc = leaf.name
            entry = tax_mapping.loc[acc]
            leaf.add_feature('tax_superkingdom', entry['superkingdom'])

            leaf.add_feature('tax_class', entry['class'])

            if leaf.tax_superkingdom == 'Eukaryota':
                leaf.add_feature('tax_filter', 'Eukaryota')
            else:
                leaf.add_feature('tax_filter', entry['class'])

        except KeyError:
            leaf.add_feature('superkingom', 'ERROR')
            leaf.add_feature('class', 'ERROR')


# return closest non self leaf
def get_closest_leaf(leaf):
    near_leaves = [near_leaf for near_leaf in leaf.up.get_leaves() if near_leaf != leaf]
    distances = [leaf.get_distance(near_leaf) for near_leaf in near_leaves]
    min_dist = min(distances)
    closest_leaf = near_leaves[distances.index(min_dist)]

    return closest_leaf, min_dist




# starting from one leaf with an attribute traverse upwards untill
# all leaves from the ancestor is no loger monophyletic under the given attribute
# repeat for all remaining leaves
# if any clade would have the global root as ancestor rerooot and retry to avoid false paraphyly by tree data struture
def get_paraphyletic_groups(tree, attribute='tax_superkingdom', attr_value='Eukaryota', current_root=False):
    # tree.set_outgroup(tree.get_farthest_leaf()[0])

    if current_root:
        tree.set_outgroup(current_root)
    else:
        current_root = tree.get_tree_root()

    # get a list of all leaves with an attribute matching the match value provided
    check_leaves = [leaf for leaf in tree.get_leaves() if getattr(leaf, attribute) == attr_value]
    clade_nodes = []

    seed_node = check_leaves[0]

    while check_leaves:
        # assume monophyly
        mono = True

        # check for all parent leaves if attribute matches the value, if not its not monophyletic, break
        for leaf in seed_node.up.get_leaves():
            if getattr(leaf, attribute) != attr_value:
                mono = False
                break

        # if monophyletic try higher node
        if mono:
            seed_node = seed_node.up

        # else retrun node and exclude all leaves from list of leaves to check
        else:
            clade_nodes.append(seed_node)
            check_leaves = [leaf for leaf in check_leaves if leaf not in seed_node.get_leaves()]
            if check_leaves:
                seed_node = check_leaves[0]

                # if parent has no parent it is the root
    if [node for node in clade_nodes if node.up.up == None]:
        # print('A tree clade has rooted parent nodes, rerooting')

        # get the first non-clade daughter from current root
        non_clade_daughter = [node for node in current_root.children if node not in clade_nodes][0]

        return get_paraphyletic_groups(tree, attribute=attribute, attr_value=attr_value,
                                       current_root=non_clade_daughter)

    else:
        return clade_nodes


# return the entropy of decendant and non decendant leaf labels
def get_entropy_for_partition(tree, node, attribute='tax_filter', attr_value='Eukaryota'):
    all_labels = [getattr(leaf, attribute) for leaf in tree.get_leaves()]
    all_label_count = all_labels.count(attr_value)
    tree_width = len(all_labels)

    #     base_label_Px = (all_label_count/tree_width)
    #     base_label_H = base_label_Px*np.log2(1/base_label_Px)

    clade_labels = [getattr(leaf, attribute) for leaf in node.get_leaves()]
    clade_label_count = clade_labels.count(attr_value)
    clade_width = len(clade_labels)

    # calculate label entropy
    # print('AAA', clade_labels, clade_label_count, clade_width)
    label_Px = (clade_label_count) / (clade_width)
    label_H = label_Px * np.log2(1 / label_Px)

    # calculate external entropy change

    # if all labels in the clade the external entropy is 0
    if all_label_count - clade_label_count == 0:
        external_label_H = 0

    else:
        # calculate external entropy change
        external_label_Px = (all_label_count - clade_label_count) / (tree_width - clade_width)
        external_label_H = external_label_Px * np.log2(1 / external_label_Px)

    return label_H, external_label_H


# assign soft LCA node based on minimizing entropy between given label outside and inside clade
# more pessimissive than voting ratio, qualitatively underestimates
def get_soft_LCA_by_relative_entropy(tree, attribute='tax_superkingdom', attr_value='Eukaryota', save_loss=False):
    # count vote ratio for each node for one taxa
    tree_width = len(tree.get_leaves())
    root = tree.get_tree_root()

    lowest_total_H = float("inf")
    best_node = root
    vote_label = attr_value

    all_labels = [getattr(leaf, attribute) for leaf in tree.get_leaves()]
    all_label_count = all_labels.count(attr_value)
    tree_width = len(all_labels)

    all_label_Px = (all_label_count / tree_width)
    all_labels_H = all_label_Px * np.log2(1 / all_label_Px)

    # for debugging
    if save_loss:
        for node in tree.traverse():
            node.add_feature('vote_loss', 'None')

    # check for monophyly
    LCA_groups = get_paraphyletic_groups(tree, attribute=attribute, attr_value=attr_value)
    if len(LCA_groups) == 1:
        print(f'The attribute {attribute} is monophyletic for {attr_value}. Returning LCA node.')
        return (LCA_groups[0], 0)

    # the best partition will be on the path from an LCA node to the root
    tested_nodes = []
    for node in LCA_groups:

        node_label_count = len([leaf for leaf in node.get_leaves() if getattr(leaf, attribute) == attr_value])

        # ascend until all labeled leaves are decendants of node
        while node != root:

            # skip known nodes
            if node not in tested_nodes:

                # calculate internal and external entropy
                label_H, external_label_H = get_entropy_for_partition(tree, node, attribute=attribute,
                                                                      attr_value=attr_value)
                total_H = label_H + external_label_H

                # penalize leaf LCAs to avoid laddered LCAs when having repeated outgroups
                # not neccesarily waned as spread singleons get their global LCA as soft_LCA
                if node.is_leaf():
                    total_H += 0.5

                # update best guess
                if total_H < lowest_total_H:
                    lowest_total_H = total_H
                    best_node = node

                if save_loss:
                    node.add_feature('vote_loss', total_H)

            # break after calculations if all nodes are decendants
            node_label_count = len([leaf for leaf in node.get_leaves() if getattr(leaf, attribute) == attr_value])
            if node_label_count == all_label_count:
                break
            # print(node_label_count, all_label_count)

            # ascend
            tested_nodes.append(node)
            node = node.up

    # print(f'Best node for {attr_value} has a total H of {lowest_total_H}')
    return (best_node, lowest_total_H)


# find soft LCA, mask all members and repeat LCA check until no more members
# returns list of all soft LCAs
def get_multiple_soft_LCAs(tree, attribute='tax_filter', attr_value='Eukaryota', min_size=1, min_purity=0,
                           max_entropy=9999):
    print(f'Searching for LCA_nodes by checking where {attribute} is {attr_value}')

    soft_LCA_nodes = []
    total_nodes = len(
        [getattr(node, attribute) for node in tree.get_leaves() if getattr(node, attribute) == attr_value])

    while total_nodes > 0:

        soft_LCA_node, lowest_H_loss = get_soft_LCA_by_relative_entropy(tree, attribute=attribute,
                                                                        attr_value=attr_value)

        soft_LCA_node_leaves = soft_LCA_node.get_leaves()
        soft_LCA_node_size = len(soft_LCA_node_leaves)

        # mark labeled included nodes and decrement total nodes left to check
        for node in soft_LCA_node_leaves:
            if getattr(node, attribute) == attr_value:
                setattr(node, attribute, 'SAMPLED')
                total_nodes -= 1

        soft_LCA_nodes.append(soft_LCA_node)

    # reset leaf node attributes
    for node in tree.get_leaves():
        if getattr(node, attribute) == 'SAMPLED':
            setattr(node, attribute, attr_value)

    # recalculate entropies
    print(f'\tRecalculating sizes, purities and entropies for all LCA nodes')
    filtered_soft_LCA_nodes = []

    for i, node in enumerate(soft_LCA_nodes):
        label_H, external_label_H = get_entropy_for_partition(tree, node, attribute=attribute, attr_value=attr_value)

        soft_LCA_members = [getattr(leaf, attribute) for leaf in node.get_leaves()]
        soft_LCA_size = len(soft_LCA_members)
        soft_LCA_purity = soft_LCA_members.count(attr_value) / soft_LCA_size
        soft_LCA_entropy = label_H + external_label_H

        if soft_LCA_size >= min_size and soft_LCA_purity >= min_purity and soft_LCA_entropy <= max_entropy:
            print(
                f'\tFound node of size {soft_LCA_size} with purity of {soft_LCA_purity} and entropy {soft_LCA_entropy} as LCA for {attr_value}')
            filtered_soft_LCA_nodes.append([node, soft_LCA_size, soft_LCA_purity, soft_LCA_entropy])
        else:
            print(
                f'\tRejected node of size {soft_LCA_size} with purity of {soft_LCA_purity} and entropy {soft_LCA_entropy} as LCA for {attr_value}')

    if len(filtered_soft_LCA_nodes) == 0:
        print(
            f'WARNING: No valid LCA nodes present for {attribute} = {attr_value} under conditions min_size={min_size}, min_purity={min_purity}, max_entropy={max_entropy}')

    print()

    return filtered_soft_LCA_nodes
