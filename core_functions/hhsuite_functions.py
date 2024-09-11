#format directory containing ONLY .fasta MSAs into hhsuite database
def build_HHsuite_database_from_MSAs(fasta_root, output_root, DB_name, mpis):

    import os
    import subprocess
    from paths_and_parameters import hhconsensus_opts, hhmake_opts, cstranslate_opts

    os.environ['OMP_NUM_THREADS'] = '1'

    output_name = f'{output_root}{DB_name}'
    tmp_name = f'{output_name}_tmp'

    mpi_setup = f'mpirun -np {mpis}'
    hhconsensus_command = f'hhconsensus -i stdin -oa3m stdout {hhconsensus_opts}'
    hhmake_command = f'hhmake -i stdin -o stdout {hhmake_opts}'

    run_dict = {
        #build sorted ffindex
        'run_ffindex' : f'ffindex_build -s {tmp_name}.ffdata {tmp_name}.ffindex {fasta_root}',

        #make a3m alignments and add consensus sequences
        'run_hhconsensus' : f'{mpi_setup} ffindex_apply_mpi {tmp_name}.ffdata {tmp_name}.ffindex -i {tmp_name}_a3m.ffindex -d {tmp_name}_a3m.ffdata -- {hhconsensus_command}',

        #create hhm profiles
        'run_hhmake' : f'{mpi_setup} ffindex_apply_mpi {tmp_name}_a3m.ffdata {tmp_name}_a3m.ffindex -i {tmp_name}_hhm.ffindex -d {tmp_name}_hhm.ffdata -- {hhmake_command}',

        #add context states
        'run_cstranslate' : f'{mpi_setup} cstranslate_mpi {cstranslate_opts} -i {tmp_name}_a3m -o {tmp_name}_cs219',

        #reorder and rename all output databases
        'run_calculate_sort' : f'sort -k3 -n -r {tmp_name}_cs219.ffindex | cut -f1 > {tmp_name}.sort',
        'run_apply_sort_hhm' : f'ffindex_order {tmp_name}.sort {tmp_name}_hhm.ffdata {tmp_name}_hhm.ffindex {tmp_name}_hhm_ordered.ffdata {tmp_name}_hhm_ordered.ffindex',
        'run_apply_sort_a3m' : f'ffindex_order {tmp_name}.sort {tmp_name}_a3m.ffdata {tmp_name}_a3m.ffindex {tmp_name}_a3m_ordered.ffdata {tmp_name}_a3m_ordered.ffindex',
        'run_rename_a3m_data' : f'mv {tmp_name}_a3m_ordered.ffindex {output_name}_a3m.ffindex',
        'run_rename_a3m_index' : f'mv {tmp_name}_a3m_ordered.ffdata {output_name}_a3m.ffdata',
        'run_rename_hhm_data' : f'mv {tmp_name}_hhm_ordered.ffindex {output_name}_hhm.ffindex',
        'run_rename_hhm_index' : f'mv {tmp_name}_hhm_ordered.ffdata {output_name}_hhm.ffdata',
        'run_rename_cs219_data': f'mv {tmp_name}_cs219.ffindex {output_name}_cs219.ffindex',
        'run_rename_cs219_index': f'mv {tmp_name}_cs219.ffdata {output_name}_cs219.ffdata',

        #clean
        'run_clean' : f'rm {tmp_name}*'
    }

    for key, value in run_dict.items():
        print(key, value, sep='\n')
        subprocess.run(value, shell=True)

    return output_name

#format and run hhsuite search swarm query and targetDB is root name as XXX (without _hhm.ffdata)
def hhsuite_swarm_search(queryDB, targetDB, output_root, splits, swarm_opts, threads=1, search_type='hhblits', search_opts='', run_with_lscratch=True):

    import subprocess

    from paths_and_parameters import path_tmp, exe_hhsuite_search, exe_python
    from core_functions.helper_functions import swarm_submit_with_lock

    #create output folders
    path_splits = queryDB+'_tmp_splits/'
    basename = queryDB.split('/')[-1]
    subprocess.run(f'mkdir -p {path_splits}'.split())


    # format swarm file
    swarmfile = path_splits + 'hhsuite_search.swarm'
    swarm = open(swarmfile, 'w')

    # format swarm header for submission
    swarm_opts['threads-per-process'] = threads
    swarm.write(''.join([f'#SWARM --{key} {value}\n' for key, value in swarm_opts.items()]))


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
        swarm.write(f'{exe_python} -u {exe_hhsuite_search} --query_hhm_ffindex {index_batch_name} --query_hhm_ffdata {queryDB}_hhm.ffdata --targetDB {targetDB} --output_root {output_root} --threads {threads} --search_type {search_type} --search_opts "{search_opts}" --run_with_lscratch {run_with_lscratch} \n')

    swarm.close()

    swarm_submit_with_lock(swarmfile, refresh=60)

    #remove temporary files
    subprocess.run(f'rm -r {path_splits}'.split())

    return output_root

#OLD
def parse_hhsuite_search_results(search_ffdata):

    import core_functions.parseHHsuite as HH

    basename = search_ffdata.split('.ffdata')[0]
    output_pkl = basename+'.pkl'
    output_pkl_filtered = basename+'.filtered.pkl'
    output_tsv_filtered = basename + '.filtered.tsv'

    #parse ffdata and save .pkl
    new_data = HH.load_HHBlitsData(search_ffdata)
    new_data.write_pkl(output_pkl)

    print(f'Parsing as HHBlitsData object: {search_ffdata}')
    for key, query in new_data.data.items():
        query.add_self_hit()
        #query.filter_numeric(field='Pairwise_cov', min=20, replace=True, keep_self=True)
        query.filter_numeric(field='Prob', min=50, replace=True, keep_self=True)

    print(f'Writing {output_pkl_filtered}')
    new_data.write_pkl(output_pkl_filtered)
    new_data.write_data_tsv(output_tsv_filtered)

    return

#parse and merge all search.ffdata from a folder to a single parse_HHsuite object
#FILTERS ARE BROKEN!?
def merge_hhsuite_search_results(search_root, output_basename, write_tsv=True, filter_cov=False, filter_prob=False):

    import os
    import core_functions.parseHHsuite as HH

    print(f'Merging all files from {search_root} to {output_basename}')

    #find all files with .ffdata
    ffdata_files = [search_root+file for file in os.listdir(search_root) if file.endswith('ffdata')]

    #initialise output object
    all_data = HH.HHblitsData()

    #loop over files and add to output
    for i, file in enumerate(ffdata_files):
        print(f'Parsing file #{i}: {file}')
        #load new data object from file
        new_data = HH.load_HHBlitsData(file)

        #filter new data object
        if filter_cov or filter_prob:
            for key, query in new_data.data.items():
                query.add_self_hit()

                if filter_cov:
                    query.filter_numeric(field='Pairwise_cov', min=filter_cov, replace=True, keep_self=True)

                if filter_prob:
                    query.filter_numeric(field='Prob', min=filter_prob, replace=True, keep_self=True)

        #merge with all_data
        all_data.add_entries(new_data.data)

    #all_data.write_pkl(output_basename+'.pkl')

    if write_tsv:
        all_data.write_query_tsv(output_basename+'.query.tsv')
        all_data.write_data_tsv(output_basename+'.tsv')
        os.system(f"rm -rf {search_root}")

    return all_data

#parse and merge all search.blast.ffdata from a folder to a single tsv
def merge_hhsuite_blasttab_results(search_root, output_basename):

    import os

    #find all files with .ffdata
    ffdata_files = [search_root+file for file in os.listdir(search_root) if file.endswith('blast.ffdata')]

    with open(output_basename+'.tsv', 'a') as outfile:

        #loop over files and add to output
        for i, file in enumerate(ffdata_files):
            print(f'Parsing file #{i}: {file}')
            with open(file, 'r') as infile:
                outfile.write(infile.read().replace('\x00', ''))

    return output_basename+'.tsv'




















