
import argparse


#prepare all mmseqs
def microcosm_run(file_root, file_basename, threads, save_intermediate_files):

    import subprocess
    import pandas as pd
    import core_functions.microcosm_functions as microcosm
    from paths_and_parameters import microcosm_format_opts
    from core_functions.software_wrappers import muscle_ensamble, calculate_IQtree
    from core_functions.helper_functions import fasta_to_dict, filter_by_entropy, dict_to_fasta
    from core_functions.tree_functions import format_leafDF

    print(f'STARTED analysing {file_root}')
    if save_intermediate_files:
        print('Keeping intermediate files')
    else:
        print('DELETING intermediate files')

    if threads == -1:
        threads = microcosm_format_opts['threads']

    max_prok_sequences = microcosm_format_opts['max_prok_sequences']
    max_euk_sequences = microcosm_format_opts['max_euk_sequences']

    print(f'Reading taxonomy mapping from {microcosm_format_opts["taxonomy_mapping"]}')


    #scrub DELETE entries from microcosms
    taxonomy_mapping = pd.read_csv(microcosm_format_opts['taxonomy_mapping'], sep='\t',
                                   names=['acc', 'orgid', 'superkingdom', 'class'], index_col=0)

    # query accs
    accs = pd.read_csv(file_root + file_basename + '.query.acc', header=None, names=['acc'], index_col=0)
    filter_accs = taxonomy_mapping.loc[accs.index] != 'DELETE'
    index_to_keep = filter_accs[filter_accs['class']].index
    pd.Series(index_to_keep).to_csv(file_root + file_basename + '.query.acc', index=None, header=None)

    # target accs
    accs = pd.read_csv(file_root + file_basename + '.target.acc', header=None, names=['acc'], index_col=0)
    filter_accs = taxonomy_mapping.loc[accs.index] != 'DELETE'
    index_to_keep = filter_accs[filter_accs['class']].index
    pd.Series(index_to_keep).to_csv(file_root + file_basename + '.target.acc', index=None, header=None)

    #preapare fasta files from accession lists
    query_fasta, target_fasta = microcosm.prepare_mmseqs(file_root,
                                                         file_basename,
                                                         original_query_DB=microcosm_format_opts['original_query_DB'],
                                                         original_target_DB=microcosm_format_opts['original_target_DB'],
                                                         save_intermediate_files=save_intermediate_files)

    #reduce both query and target to size using FAMSA and Fasttree

    # if there are too many query sequences, crop to size
    with open(query_fasta, 'r') as fastafile:
        size = fastafile.read().count('>')

    if size > max_euk_sequences:
        print(f"There are more than {max_euk_sequences} sequences in {query_fasta} ({size}), will crop considering taxonomy in {microcosm_format_opts['taxonomy_mapping']}")

        query_fasta = microcosm.fasta_reduce_size(query_fasta,
                                                  threads=threads,
                                                  max_leaf_size=max_euk_sequences,
                                                  filter_entropy=microcosm_format_opts['filter_entropy'],
                                                  save_intermediate_files=True,
                                                  taxDF=taxonomy_mapping,
                                                  predelete_outliers=True)

    else:
        print(f' There are less than {max_euk_sequences} sequences in {query_fasta} ({size}), no cropping needed')
        print(f' Writing direct .leaf_mapping file for {query_fasta} as no cropping is done.')

        with open(file_root + file_basename + '.query.acc', 'r') as acc_file:
            accs = [acc.strip() for acc in acc_file.readlines()]

        crop_dict = {acc:[acc] for acc in accs}
        leafDF = format_leafDF(crop_dict, taxonomy_mapping)
        leafDF.to_csv(query_fasta + '.leaf_mapping', header=None, sep='\t')


    # if there are too many target sequences, crop to size
    with open(target_fasta, 'r') as fastafile:
        size = fastafile.read().count('>')

    if size > max_prok_sequences:
        print(f'There are more than {max_prok_sequences} sequences in {target_fasta} ({size}), will crop considering taxonomy in {microcosm_format_opts["taxonomy_mapping"]}')

        target_fasta = microcosm.fasta_reduce_size(target_fasta,
                                                   threads = threads,
                                                   max_leaf_size=max_prok_sequences,
                                                   filter_entropy = microcosm_format_opts['filter_entropy'],
                                                   save_intermediate_files=True,
                                                   taxDF=taxonomy_mapping,
                                                   predelete_outliers=True)

    else:
        print(f' There are less than {max_prok_sequences} sequences in {target_fasta} ({size}), no cropping needed')
        print(f' Writing direct .leaf_mapping file for {target_fasta} as no cropping is done.')

        with open(file_root + file_basename + '.target.acc', 'r') as acc_file:
            accs = [acc.strip() for acc in acc_file.readlines()]

        crop_dict = {acc:[acc] for acc in accs}
        leafDF = format_leafDF(crop_dict, taxonomy_mapping)
        leafDF.to_csv(target_fasta + '.leaf_mapping', header=None, sep='\t')


    # merge fasta files and leaf mapping and realign
    merged_fasta = file_root + file_basename + '.merged.fasta'
    subprocess.run(f'cat {query_fasta} {target_fasta} > {merged_fasta}', shell=True)
    #subprocess.run(f'cat {file_root + file_basename}  {file_root + file_basename} > {merged_fasta}', shell=True)
    subprocess.run(f'cat {query_fasta + ".leaf_mapping"} {target_fasta + ".leaf_mapping"} > {merged_fasta + ".leaf_mapping"}', shell=True)

    print(f'Merged prok and euk files and aligning using muscle')
    muscle_fasta = muscle_ensamble(merged_fasta,
                    threads = threads,
                    muscle_reps = microcosm_format_opts['muscle_reps'],
                    super5=False,
                    save_intermediate_files=save_intermediate_files,
                    timeout=microcosm_format_opts['muscle_timeout'])

    # filter by entropy
    aln = fasta_to_dict(file=muscle_fasta)
    aln_filter = filter_by_entropy(seq_dict = aln,
                                   entropy_min = microcosm_format_opts['filter_entropy'])

    if save_intermediate_files:
        subprocess.run(f'cp {muscle_fasta} {muscle_fasta}.uncropped'.split())

    dict_to_fasta(aln_filter, write_file=muscle_fasta)

    #create IQTree from alignment
    treefile = calculate_IQtree(muscle_fasta,
                                evo_model_params = microcosm_format_opts['evo_model_params'],
                                threads = threads,
                                save_intermediate_files=save_intermediate_files)

    # process final tree, find outliers, find LCAs for all taxa and calculate distances
    annot_tree, tree_data = microcosm.tree_analysis(treefile,
                                              merged_fasta + ".leaf_mapping",
                                              file_basename,
                                              delete_outliers=True)

    # write to file before
    tree_data.to_csv(file_root + file_basename + '.merged.tree_data.tsv', sep='\t')

    # draw final master tree and write to file
    #color_tree = microcosm.color_tree(annot_tree)
    annot_tree.write(features=['LCA', 'counts', 'dist', 'name', 'support', 'taxa'], outfile=treefile + '.annot',
                     format=1)

    # break in cases of poor EUK annotation
    if tree_data.shape[0] == 0:
        print(f'No tree data was detected in {file_basename} exiting')
        return

    # skip constraint analysis if only one prok taxa found
    if tree_data.prok_clade_rep.unique().shape[0] == 1:

        print(f'Only one prok LCA was detected, exiting without clade analysis')

        # add fake lines to indicate absolute certainty hit
        test_data_cols = ['logL', 'deltaL', 'bp-RELL', 'bp-RELL_accept', 'p-KH', 'p-KH_accept',
                          'p-SH', 'p-SH_accept', 'c-ELW', 'c-ELW_accept', 'p-AU', 'p-AU_accept']

        empty_test_data = pd.DataFrame('+', index=range(tree_data.shape[0]), columns=test_data_cols)
        empty_test_data[['bp-RELL', 'p-KH', 'p-SH', 'c-ELW', 'p-AU']] = [[1] * 5] * tree_data.shape[0]
        empty_test_data[['logL', 'deltaL']] = [[-1] * 2] * tree_data.shape[0]

        merged_data = pd.concat([tree_data.reset_index(), empty_test_data], axis=1)
        merged_data.to_csv(file_root + file_basename + '.merged.tree_data.tsv', sep='\t')

        return


    # format constraint trees for clade likelihood analysis
    constraint_job_data = microcosm.format_constraint_analysis(file_root, file_basename, tree_data)

    # dirty grep of model selected by modelfinder
    with open(file_root + file_basename + '.merged.fasta.muscle.iqtree', 'r') as iqtreefile:
        best_fit_text = 'Best-fit model according to'
        evo_model = [line.rsplit(':',1)[1].strip() for line in iqtreefile.readlines() if line.startswith(best_fit_text)][0]

    # calculate constrained IQtrees and format results, merge into tree_data
    all_test_data = microcosm.run_constraint_analysis(constraint_job_data, evo_model, threads)

    test_data_cols = ['logL', 'deltaL', 'bp-RELL', 'bp-RELL_accept', 'p-KH', 'p-KH_accept',
                     'p-SH', 'p-SH_accept', 'c-ELW','c-ELW_accept', 'p-AU', 'p-AU_accept']

    #reset indexes and
    merged_data = pd.concat([tree_data.sort_values(by=['euk_clade_rep', 'prok_clade_rep']).reset_index(),
                             all_test_data.sort_values(by=['euk_clade_rep', 'prok_clade_rep']).reset_index().loc[:,test_data_cols]], axis=1)
    merged_data.to_csv(file_root + file_basename + '.merged.tree_data.tsv', sep='\t')

    # do cleanup
    if save_intermediate_files == False:
        subprocess.run(f'rm {query_fasta} {target_fasta}'.split())
        subprocess.run(f'rm {merged_fasta}.muscle-efa {merged_fasta}.muscle.log'.split())

        subprocess.run(f'rm {file_root}/*.fasttree {file_root}/*muscle.log {file_root}/*.famsa* {file_root}/*.uncropped', shell=True)
        subprocess.run(f'rm {file_root}/*.contree {file_root}/*.bionj {file_root}/*.gz {file_root}/*.mldist {file_root}/*.nex {file_root}/*.ufboot', shell=True)


#argparse define
parser = argparse.ArgumentParser(description='Evaluate microcosm')
parser.add_argument('--file_basename', type=str, required=True, help='query root name, maps to folder containing containing .acc and .targets accession files within microcosm')
parser.add_argument('--file_root', type=str, required=True, help='path to microcosm root')
parser.add_argument('--threads', type=int, nargs='?', default=-1, const=-1, help='path to microcosm root')
parser.add_argument('--delete_intermediate_files', action='store_false', help='remove all intermediate files')

args = parser.parse_args()

#run main
if __name__ == '__main__':

    microcosm_run(args.file_root,
                  args.file_basename,
                  args.threads,
                  save_intermediate_files=args.delete_intermediate_files)















