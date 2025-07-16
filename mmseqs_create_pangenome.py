# +
import argparse
import subprocess
import pandas as pd

# requires mmseqs databse with preannotated taxonomy
# produced mmseqsDB with taxonomy
def mmseqs_form_pangenome(mmseqs_seqDB,
                   pangenome_seqBD,
                   evaluation_rank = 'family',
                   evaluation_cutoff = 0.5,
                   threads=16,
                   mmseqs_exe = 'mmseqs',
                   mmseqs_cluster_opts = '-s 3 -c 0.8 --cov-mode 0',
                   save_intermedidate_files = True,
                   ):
    
    # housekeeping and input checking
    fam_abbr = {'superkingdom':'d', 'phylum':'p', 'class':'c', 'order':'o',
                'family':'f', 'genus':'g', 'species':'s'}
    
    if evaluation_rank not in fam_abbr.keys():
        print(f'Warning! Evaluation_rank: {evaluation_rank} not in {fam_abbr.keys()}')
        return

    #define paths
    working_root  = pangenome_seqBD + '_tmp'
    threads = f'--threads {threads}'
    
    # make working directory structure
    subprocess.run(f'mkdir {working_root}'.split())
    subprocess.run(f'mkdir {working_root}/mmseqs_tmp'.split())

    # msmeqs calculations
    # calculate taxonomy
    subprocess.run(f'{mmseqs_exe} taxonomyreport {mmseqs_seqDB} {mmseqs_seqDB} {working_root}/taxreport '.split())
    
    # run cluster for each DB agaisnt itself -s 7 -e 1E-10 -c 0.8 --cov-mode 0 no min seq id
    subprocess.run(f'{mmseqs_exe} cluster {mmseqs_seqDB} {working_root}/cluster {working_root}/mmseqs_tmp {mmseqs_cluster_opts} {threads}'.split())
    
    # add taxonomy with --tax-lineage 1 to get lineages
    subprocess.run(f'{mmseqs_exe} addtaxonomy --tax-lineage 1 {mmseqs_seqDB} {working_root}/cluster {working_root}/cluster {threads}'.split())
    
    # make tsv for parsing
    subprocess.run(f'{mmseqs_exe} createtsv {mmseqs_seqDB} {mmseqs_seqDB} {working_root}/cluster {working_root}/cluster.tsv'.split())


    # results parsing
    # load tax report to find number of tax denominations
    taxDF = pd.read_csv(f'{working_root}/taxreport',
                                  sep='\t', header=None, index_col=4,
                                  names = ['frac', 'sum_counts','counts','rank','taxid','taxon'],
                                  dtype={'frac':float, 'sum_counts':int,'counts':int,'rank':str,'taxid':int,'taxon':str})
    taxDF.taxon = taxDF.taxon.str.strip()
    
    # calculate number of valid taxa for given rank
    all_valid_ranks = taxDF[taxDF['rank'] == evaluation_rank].taxon.values
    num_valid_ranks = len(set(all_valid_ranks))

    
    # load search data
    cluster_data = pd.read_csv(f'{working_root}/cluster.tsv', sep = '\t', 
                              names = ['query','target','taxid','rank','taxname','lineage'], index_col=0)
    
    whitelist = []
    select_clusters = []
    unique_valid_ranks = []
    
    for i, data in cluster_data.groupby('query', sort=False):
        
        all_existing_ranks = set(entry[2:] for lin in data.lineage for entry in lin.split(';') if entry[0] == fam_abbr[evaluation_rank])

        unique_valid_ranks.extend([len(all_existing_ranks)]*data.shape[0])

    # add tax count and save to .tsv
    cluster_data['unique_valid_ranks'] = unique_valid_ranks
    cluster_data.reset_index().to_csv(f'{working_root}/cluster.tsv', sep='\t', index=None)
    
    # filter and save data which meet inclusion criteria
    valid_data = cluster_data[cluster_data.unique_valid_ranks > num_valid_ranks*evaluation_cutoff]['target']
    valid_data.to_csv(f'{working_root}/valid_accs', sep='\t', index=None, header=None)
    
    # create new subdb from accessions
    subprocess.run(f'{mmseqs_exe} createsubdb --subdb-mode 0 --id-mode 1 {working_root}/valid_accs {mmseqs_seqDB} {pangenome_seqBD}'.split())

    if not save_intermedidate_files:
        subprocess.run(f'rm -r {working_root}'.split())
    
    return pangenome_seqBD
    

#argparse define
parser = argparse.ArgumentParser(description='Cluster, evaluate taxonomic coverage and extract pangenome using mmseqs',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--seqDB', type=str, required=True, help='mmseqsDB with taxonomy information')
parser.add_argument('--outputDB', type=str, required=True, help='resulting pangenome mmseqsDB')
parser.add_argument('--threads', type=int, required=True, help='threads to run for all subprocesses')
parser.add_argument('--mmseqs_cluster_opts',type=str,  required=False,  nargs='?', const=1, default='-s 3 -c 0.8 --cov-mode 0', help='optional parameters to pass to mmseqs cluster for pangenome cluster creation')
parser.add_argument('--evaluation_rank', type=str, required=True, help='NCBI taxonomic tank from which to evaluate pangenome clusters over. Choose from "superkingdom, phylum, class, order, family, genus, species')
parser.add_argument('--evaluation_cutoff', type=float, required=True, help='fraction of total ranks which must be present within pangenome clusters')
parser.add_argument('--mmseqs_exe', type=str, nargs='?', const=1, default='mmseqs', help='mmseqs executable path')
parser.add_argument('--save_intermediate_files', type=bool, nargs='?', const=1, default=True, required=False)
args = parser.parse_args()

#run main
if __name__ == '__main__':
    
    mmseqs_form_pangenome(mmseqs_seqDB = args.seqDB, 
                       pangenome_seqBD = args.outputDB,
                       evaluation_rank = args.evaluation_rank, 
                       evaluation_cutoff = args.evaluation_cutoff, 
                       threads = args.threads, 
                       mmseqs_exe = args.mmseqs_exe, 
                       mmseqs_cluster_opts = args.mmseqs_cluster_opts,
                       save_intermedidate_files = args.save_intermediate_files
                       )


