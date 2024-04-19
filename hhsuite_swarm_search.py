import os
import argparse
from uuid import uuid4
import subprocess

from paths_and_parameters import exe_hhblits_omp, exe_hhsearch_omp

#format and run hhsuite search swarm query and targetDB is root name as XXX (without _hhm.ffdata)
def hhsuite_swarm_search(queryDB, targetDB, output_root, splits, threads, search_type='hhblits', run_with_lscratch=True):
    import subprocess

    from paths_and_parameters import path_tmp, exe_hhsuite_search, exe_python, hhblits_swarm_opts, hhblits_search_opts
    from core_functions.helper_functions import swarm_submit_with_lock

    #create output folders
    path_splits = queryDB+'_tmp_splits/'
    basename = queryDB.split('/')[-1]
    subprocess.run(f'mkdir -p {path_splits}'.split())


    # format swarm file
    swarmfile = path_splits + 'hhsuite_search.swarm'
    swarm = open(swarmfile, 'w')

    # format swarm header for submission
    hhblits_swarm_opts['threads-per-process'] = threads
    swarm.write(''.join([f'#SWARM --{key} {value}\n' for key, value in hhblits_swarm_opts.items()]))


    #read indexfile for query
    with open(queryDB+'_hhm.ffindex', 'r') as queryfile:
        queries = queryfile.readlines()

    # python ceiling division
    split_size = -(len(queries)//-splits)

    # partition query index into batches of max_size, write to temporary dir
    for split in range(0, len(queries), max(split_size, 1)):
        index_batch = queries[split:split + split_size]
        index_batch_name = f'{path_splits}{basename}_{split//split_size}_hhm.ffindex'

        #write index batch
        with open(index_batch_name, 'w') as batchfile:
            batchfile.writelines(index_batch)

        #add execution line to swarmfile
        swarm.write(f'{exe_python} -u {exe_hhsuite_search} --query_hhm_ffindex {index_batch_name} --query_hhm_ffdata {queryDB}_hhm.ffdata --targetDB {targetDB} --output_root {output_root} --threads {threads} --search_type {search_type} --search_opts "{hhblits_search_opts}" --run_with_lscratch {run_with_lscratch} \n')

    swarm.close()

    swarm_submit_with_lock(swarmfile, refresh=60)

    #remove temporary files
    subprocess.run(f'rm -r {path_splits}'.split())

    return output_root

#argparse define
parser = argparse.ArgumentParser(description='search query hhsuite (style) db against target hhsuite db, splitting queries and searching in parallel on lscratch copies of the target')
parser.add_argument('--queryDB', type=str, required=True, help='input full query database path name </data/user/.../XXX>_hhm.{ffindex,ffdata}')
parser.add_argument('--targetDB', type=str, required=True, help='input full target database path name </data/user/.../YYY>hhm.{ffindex,ffdata}')
parser.add_argument('--output_root', type=str, required=True, help='hhs file output root folder')
parser.add_argument('--splits', type=int, required=True, help='number of swarm jobs')
parser.add_argument('--threads', type=int, required=False, nargs='?', default=1, help='OMP threads to run')
parser.add_argument('--search_type', type=str, required=False, nargs='?', default='hhblits', help='select from "hhsearch" or "hhblits"')
parser.add_argument('--run_with_lscratch', type=bool, required=False, nargs='?', default=True, help='this tag move temporary files to biowulf /lscratch/$SLURM_JOB_ID/hhsuite_search.uuid')
args = parser.parse_args()


#run main
if __name__ == '__main__':
    hhsuite_swarm_search(queryDB=args.queryDB,
                   targetDB = args.targetDB,
                   output_root = args.output_root,
                   splits = args.splits,
                   threads = args.threads,
                   search_type = args.search_type,
                   run_with_lscratch=args.run_with_lscratch
                   )