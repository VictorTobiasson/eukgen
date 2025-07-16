# +
# This should be extensively refactored as the functions and file handling are really messy right now. 
# Effectively its just one protocol which does not need to be split in two functions. 
import argparse
import subprocess
import pandas as pd

# perform the partitioning and crop and realign each subcluster
def anneal_and_realign_subclusters(fasta,
                                   threads=16,
                                   max_sequences=250,
                                   min_sequences=10,
                                   filter_entropy=0,
                                   muscle_reps=5,
                                   muscle_timeout=7200,
                                   save_intermediate_files=True):
    

    
    from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy
    from core_functions.microcosm_functions import fasta_reduce_size
    from core_functions.software_wrappers import muscle_ensamble
    from Bio import SeqIO

    # since there is no root reference handle direct file references as well    
    if '/' in fasta:
        working_root = fasta.rsplit('/',1)[0]+'/'
    else:
        working_root = './'
        
    seqs = {seq.id:seq for seq in SeqIO.parse(fasta, format='fasta')}
        
    # this logic is duplicated in anneal_fasta_by_UMAP_HDBSCAN since code is not refactored
    # exit for small trees
    if len(seqs) <= min_sequences:
        print(f'No annealing as tree only has {len(seqs)} leaves')
        
        # write naive cluster mapping and continue
        leaf_names = [k for k in seqs.keys()]
        leafDF = pd.DataFrame({'acc':leaf_names, 'cluster_acc':[leaf_names[0]]*len(leaf_names)})
        leafDF.set_index('cluster_acc').to_csv(f'{fasta}.cluster.tsv', header=None, sep='\t')
        
        output_fastas = [fasta] 
        
    else:
        output_fastas, leafDF, cropped_tree = anneal_fasta_by_UMAP_HDBSCAN(fasta,                                                                        
                                                                           threads,
                                                                           max_sequences,
                                                                           min_sequences,
                                                                           filter_entropy,
                                                                           save_intermediate_files)
    
    # realign and crop resulting list of fastas
    print(f'Realigning resulting clusters')
    final_fastas = []
    for cluster_fasta in output_fastas:
        print(cluster_fasta)
        
        with open(cluster_fasta, 'r') as fastafile:
            size = fastafile.read().count('>')

        # if there are too many sequences, crop to size
        if size > max_sequences:
            print(f'There are more than {max_sequences} sequences in {cluster_fasta} ({size}), will crop to size')

            cropped_fasta = fasta_reduce_size(cluster_fasta, threads, max_sequences, filter_entropy, save_intermediate_files=save_intermediate_files)

        else:
            print(f'There are less than {max_sequences} sequences in {cluster_fasta} ({size}), no cropping needed')
            cropped_fasta = cluster_fasta

        # align
        print(f'Aligning with muscle5 as ensemble with {muscle_reps} replicates')
        muscle_fasta = muscle_ensamble(cropped_fasta, threads, muscle_reps, super5=False, 
                                       save_intermediate_files=save_intermediate_files, timeout=muscle_timeout)

        # filter by entropy
        aln = fasta_to_dict(file=muscle_fasta)
        aln_filter = filter_by_entropy(aln, filter_entropy)
        
        final_fastas.append(muscle_fasta)
        
        if save_intermediate_files:
            dict_to_fasta(aln_filter, write_file=f'{muscle_fasta}.aln')

        else:
            dict_to_fasta(aln_filter, write_file=cluster_fasta)
            subprocess.run(f'rm {cluster_fasta}.muscle {cluster_fasta}.muscle-efa {cluster_fasta}.leaf_mapping {cluster_fasta}.muscle.log'.split())
            
        
    print(f'Final aligned files are:')
    for f in final_fastas:
        print(f)
            
    return


# align and crop fasta, construct fasttree
# remove leaf and clade outlier from tree
# calculate pairwise distances from resulting tree and cluster with UMAP/HDBSCAN
def anneal_fasta_by_UMAP_HDBSCAN(merged_fasta, threads=8, max_leaves=250, min_seqs=10, filter_entropy=0, save_intermediate_files=True):

    from core_functions.tree_functions import get_outlier_nodes_by_lognorm, remove_outlier_nodes
    from core_functions.tree_functions import crop_leaves_to_size, weighted_midpoint_root, partition_tree_with_UMAP_HDBSCAN
    from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy
    from paths_and_parameters import exe_fasttree, exe_famsa

    from Bio import SeqIO
    from ete3 import Tree
    
    seqs = {seq.id:seq for seq in SeqIO.parse(merged_fasta, format='fasta')}

    print(f'Aligning and constructing FastTree for {merged_fasta} with {len(seqs)} sequences')
    

    # align, filter and construct tree for seqs
    with open(f'{merged_fasta}.famsa.log', 'a') as famsa_logfile, open(f'{merged_fasta}.fasttree.log', 'a') as fasttree_logfile:

        famsa_command = exe_famsa + f' -t {threads} {merged_fasta} {merged_fasta}.famsa'
        subprocess.run(famsa_command.split(), stdout=famsa_logfile, stderr=famsa_logfile)

        aln = fasta_to_dict(file=f'{merged_fasta}.famsa')
        aln_filter = filter_by_entropy(aln, filter_entropy)
        dict_to_fasta(aln_filter, write_file=f'{merged_fasta}.famsa')

        fasttree_command = exe_fasttree + f" -gamma -out {merged_fasta}.fasttree {merged_fasta}.famsa"
        subprocess.run(fasttree_command.split(), stdout=fasttree_logfile, stderr=fasttree_logfile)

    treefile = f'{merged_fasta}.fasttree'
    tree = Tree(treefile)
    tree = weighted_midpoint_root(tree)


    # exit for small trees
    if len(tree) <= min_seqs:
        print(f'Exiting as tree only has {len(tree)} leaves')

        # write naive cluster mapping and reassign tree reference
        leaf_names = [l.name for l in tree.get_leaves()]
        leafDF = pd.DataFrame({'acc':leaf_names, 'cluster_acc':[leaf_names[0]]*len(leaf_names)})
        leafDF.set_index('cluster_acc').to_csv(f'{merged_fasta}.cluster.tsv', header=None, sep='\t')
        
        cropped_tree = tree
     
    #otherwise crop leaves and clades, reduce tree and cluster UMAP embedding
    else:
        
        # perform operations on separate tree for safety
        cropped_tree = tree.copy()

        # remove outlier leaves and nodes to avoid retaining outliers for cropping leaves
        outlier_nodes, dist_series, fits, cutoff = get_outlier_nodes_by_lognorm(cropped_tree, p_low=0, p_high=0.90, only_leaves=True, deletion_cutoff=0.2)
        cropped_tree = remove_outlier_nodes(cropped_tree, outlier_nodes)
        print(f'After removing outliers leaves Tree has {len(cropped_tree)} leaves')

        outlier_nodes, dist_series, fits, cutoff = get_outlier_nodes_by_lognorm(cropped_tree, p_low=0, p_high=0.99, drop_leaves=True, only_leaves=False, deletion_cutoff=0.15)
        cropped_tree = remove_outlier_nodes(cropped_tree, outlier_nodes)
        print(f'After removing outliers Tree has {len(cropped_tree)} leaves')

        # crop tree to max_leaves, track remaining non-outlier leaves
        print(f'Cropping joint tree to {max_leaves}')
        cropped_tree, crop_dict = crop_leaves_to_size(cropped_tree, max_size = max_leaves, 
                                              save_cropped_references=True, monophyletic=False, crop_dict=None)

        leafDF = pd.DataFrame(pd.Series(crop_dict), columns=['acc']).explode('acc')
        leafDF.index.name = 'leaf'

        # obtain a tree partition with UMAP and HDBSCAN
        cropped_tree, UMAP_data = partition_tree_with_UMAP_HDBSCAN(cropped_tree)

        # get leaf name reference names for clusters
        set_names = {l.set:l.name for l in cropped_tree.get_leaves()}
        set_data = pd.DataFrame([[l.name, set_names[l.set], l.set] for l in cropped_tree.get_leaves()], columns = ['leaf', 'cluster_acc', 'set'])
        set_data = set_data.set_index('leaf')

        # assign all accs a cluster based on their leaf reference
        leafDF['cluster_acc'] = set_data.loc[leafDF.index].cluster_acc

        #do not save outlier sequences from HDBSCAN
        if -1 in set_names.keys():
            leafDF = leafDF[leafDF.cluster_acc != set_names[-1]]
            
    cropped_tree.write(outfile=f'{merged_fasta}.treefile.cropped.annot', features=['set'])

    # write all clusters to separate files for secondary trimming and alignment
    output_fastas = []
    print('Final cluster assignment, ignoring outliers:')
    for cluster, data in leafDF.groupby('cluster_acc'):
        accs = data.acc
        
        # since there is no root reference handle direct file references as well
        if '/' in merged_fasta:
            cropped_fasta = f'{merged_fasta.rsplit("/",1)[0]}/{cluster}'
        else:
            cropped_fasta = f'{cluster}'
            
        output_fastas.append(cropped_fasta)

        with open(cropped_fasta, 'w') as out: 
            cropped_seqs = [seqs[acc] for acc in accs]
            print(f'Cluster {cluster} has a final {len(cropped_seqs)} seqs')
            SeqIO.write(cropped_seqs, out, format='fasta')
            
    if not save_intermediate_files:
        print('Deleting intermediate files')
        subprocess.run(f'rm {merged_fasta}.famsa {merged_fasta}.famsa.log {merged_fasta}.fasttree {merged_fasta}.fasttree.log {merged_fasta}.treefile.cropped.annot'.split())
        
    # write final cluster mapping
    leafDF.set_index('cluster_acc').to_csv(f'{merged_fasta}.cluster.tsv', header=None, sep='\t')

    return output_fastas, leafDF, cropped_tree
# -

#argparse define
parser = argparse.ArgumentParser(description='Evaluate homogeneity of fasta using FastTree and cluster. Then crop subclusters to size and realign using multiple rounds of muscle')
parser.add_argument('--fasta', type=str, required=True, help='fasta file')
parser.add_argument('--threads', type=int, required=False, nargs='?', const=1, default=1, help='threads to run for all subprocesses')
parser.add_argument('--max_seqs', type=int, required=False, nargs='?', const=1, default=250, help='maximum allowed sequences without filtering')
parser.add_argument('--min_seqs', type=int, required=False, nargs='?', const=1, default=10, help='minimum cluster size to allow')
parser.add_argument('--filter_entropy', type=float, required=False, nargs='?', const=1, default=0, help='minimal column entropy to allow in alignment')
parser.add_argument('--muscle_reps', type=int, required=False, nargs='?', const=1, default=1, help='number of diversified muscle iterations, must be greater than 1')
parser.add_argument('--muscle_timeout', type=float, required=False, nargs='?', const=1, default=3600, help='soft running time timeout, if exceeded will resubmit with super5')
parser.add_argument('--save_intermediate_files', action=argparse.BooleanOptionalAction)
args = parser.parse_args()

#run main
if __name__ == '__main__':
    
    anneal_and_realign_subclusters(fasta=args.fasta,
                                   threads=args.threads,
                                   max_sequences=args.max_seqs,
                                   min_sequences=args.min_seqs,
                                   filter_entropy=args.filter_entropy,
                                   muscle_reps=args.muscle_reps,
                                   save_intermediate_files=args.save_intermediate_files,
                                   muscle_timeout=args.muscle_timeout
)
