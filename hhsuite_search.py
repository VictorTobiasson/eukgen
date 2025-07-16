import os
import argparse
from uuid import uuid4
import subprocess

from paths_and_parameters import exe_hhblits_omp, exe_hhsearch_omp

#intended to run as a split swarm submission on biowulf with lscratch as running_path
def hhsuite_search(query_hhm_ffindex, query_hhm_ffdata, targetDB, output_root, threads=1,
                   search_type='hhblits', search_opts='', run_with_lscratch=True, return_blast_m6=False):

    
    #required to avoid segfault
    os.environ['OMP_STACKSIZE'] = '32768'

    #define paths and names
    query_basename = query_hhm_ffindex.split('/')[-1].split("_hhm.ffindex")[0]
    target_basename = targetDB.split('/')[-1]
    output_name = f'{output_root}/{query_basename}_vs_{target_basename}'

    #create output root
    subprocess.run(f'mkdir -p {output_root}'.split())

    #copy query and target databse to local lscratch directory for parallell IO
    if run_with_lscratch != None:
        slurmID = os.environ['SLURM_JOB_ID']
        running_root = f'/lscratch/{slurmID}/hhsuite_search.{uuid4().hex}/'
        subprocess.run(f'mkdir -p {running_root}'.split())

    else:
        running_root = output_root

    run_dict = {
        'search_copy_data' : f'cp {query_hhm_ffdata} {running_root}{query_basename}_hhm.ffdata',
        'search_copy_index' : f'cp {query_hhm_ffindex} {running_root}',
        'search_copy_targetDB': ' '.join(['cp', targetDB+'*{ffindex,ffdata}', running_root]),
        # 'search_softlink_query_HHM_index' : f'ln -s {running_root}{query_hhm_} {running_root}{query_basename}_hhm.ffdata'
    }

    if search_type == 'hhblits':
        run_dict['search'] = f'{exe_hhblits_omp} -i {running_root}{query_basename}_hhm -d {running_root}{target_basename} -cpu {threads}, -o {output_name} {search_opts}'
        
        if return_blast_m6 != None:
            print(return_blast_m6)
            # append blast output option
            run_dict['search'] = run_dict['search'] + f' -blasttab {output_name}.blast'
        
    elif search_type == 'hhsearch':
        print('NOT IMPLEMENTED')
        #run_dict['search'] = ''
    else:
        print('"search_type" parameter must be "hhsearch" or "hhblits"')

    print(run_dict)
        
    for key, value in run_dict.items():
        print(key, value, sep='\n')
        subprocess.run(value, shell=True)

    return output_name


# +
#argparse define
parser = argparse.ArgumentParser(description='Reduce fasta to size and realign using multiple rounds of muscle')
parser.add_argument('--query_hhm_ffindex', type=str, required=True, help='input full query index name XXX_hhm.ffindex')
parser.add_argument('--query_hhm_ffdata', type=str, required=True, help='input full query data name XXX_hhm.ffindex')
parser.add_argument('--targetDB', type=str, required=True, help='input target basename "XXX", not _yyy.{ffindex,ffdata}')
parser.add_argument('--output_root', type=str, required=True, help='hhs file output root folder')
parser.add_argument('--threads', type=int, required=False, nargs='?', const=1, help='OMP threads to run')
parser.add_argument('--search_type', type=str, required=False, nargs='?', const='hhblits', help='select from "hhsearch" or "hhblits"')
parser.add_argument('--search_opts', type=str, required=False, nargs='?', const='', help='additional search parameters to pass to the search')
parser.add_argument('--run_with_lscratch', action=argparse.BooleanOptionalAction, help='move temporary files to biowulf /lscratch/$SLURM_JOB_ID/hhsuite_search.uuid')
parser.add_argument('--return_blast_m6', action=argparse.BooleanOptionalAction, help='output additional m6 format XXX.blast.{ffindex,ffdata}')

args = parser.parse_args()
# -


#run main
if __name__ == '__main__':
    
    hhsuite_search(query_hhm_ffindex=args.query_hhm_ffindex,
                   query_hhm_ffdata = args.query_hhm_ffdata,
                   targetDB = args.targetDB,
                   output_root = args.output_root,
                   threads = args.threads,
                   search_type = args.search_type,
                   search_opts = args.search_opts,
                   run_with_lscratch=args.run_with_lscratch,
                   return_blast_m6=args.return_blast_m6
                   )


