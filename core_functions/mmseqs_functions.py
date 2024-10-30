import subprocess


# extra opts is dictionary with "nonredundant_cluster": parameters
def collapse_nonredundant(root, seqDB, threads, extra_opts=None):
    from paths_and_parameters import path_mmseqs_tmp, exe_mmseqs

    # define folder structure
    subprocess.run(f'mkdir {root}/nonredundant/'.split())

    # define file names
    basename = seqDB.split('/')[-1]
    nonredundant_cluster = f'{root}nonredundant/{basename}.cluster'
    nonredundant_repseq = f'{nonredundant_cluster}.repseq'

    # define commands and parameters to run

    extra_opts_cluster = ''
    mmseqs_threads = f'--threads {threads}'

    if isinstance(extra_opts, dict):
        try:
            extra_opts_cluster = extra_opts['nonredundant_cluster']
        except KeyError:
            print(f'WARNING: Extra Opts specified for collapse_nonredundant() but no "nonredundant_cluster" found')

    # perform initial clustering for nonredundant sequences
    mmseqs_nonredundant_cluster = f'{exe_mmseqs} cluster {seqDB} {nonredundant_cluster} {path_mmseqs_tmp} {mmseqs_threads} {extra_opts_cluster}'

    # select one representative per cluster at random
    mmseqs_nonredundant_repseq = f'{exe_mmseqs} result2repseq {seqDB} {nonredundant_cluster} {nonredundant_repseq} {mmseqs_threads}'

    # c reate and return softlink for further referencing
    mmseqs_ln = f'mmseqs lndb {nonredundant_repseq} {root}/{basename}.repseq'

    # execute
    subprocess.run(mmseqs_nonredundant_cluster.split())
    subprocess.run(mmseqs_nonredundant_repseq.split())
    subprocess.run(mmseqs_ln.split())

    return f'{root}/{basename}.repseq'


# perform initial sensitive cascaded clustering, make profiles, rerun profile search;
def create_initial_profile_cluster(root, seqDB, threads, extra_opts=None):
    from paths_and_parameters import path_mmseqs_tmp, exe_mmseqs

    # define folder structure
    subprocess.run(f'mkdir {root}/initial/'.split())

    # define file names
    basename = seqDB.split('/')[-1]
    initial_cluster = f'{root}initial/{basename}.cluster'
    initial_cluster_profile = f'{initial_cluster}.profiles'
    initial_cluster_profile_consensus = f'{initial_cluster}.consensus'
    initial_search = f'{root}initial/{basename}.search'
    search_clust = f'{initial_search}.clust'
    search_clust_merged = f'{search_clust}.merged'

    # define and assign commands and parameters to run
    mmseqs_threads = f'--threads {threads}'
    extra_opts_initial_cluster = ''
    extra_opts_profile = ''
    extra_opts_search = ''
    extra_opts_clust = ''

    if isinstance(extra_opts, dict):
        try:
            extra_opts_initial_cluster = extra_opts['initial_cluster']
            extra_opts_profile = extra_opts['result2profile']
            extra_opts_search = extra_opts['initial_search']
            extra_opts_clust = extra_opts['cascade_clust']
        except KeyError:
            print(f'WARNING: Extra opts specified for initial_profile_cluster() but not all found in keys!')

    run_dict = {
        # perform first sensitive clustering and make profiles and consensus

        'mmseqs_initial_cluster': f'{exe_mmseqs} cluster {seqDB} {initial_cluster} {path_mmseqs_tmp} {mmseqs_threads} {extra_opts_initial_cluster}',
        'mmseqs_initial_cluster_r2p': f'{exe_mmseqs} result2profile {seqDB} {seqDB} {initial_cluster} {initial_cluster_profile} {mmseqs_threads} {extra_opts_profile}',
        'mmseqs_initial_p2c': f'{exe_mmseqs} profile2consensus {initial_cluster_profile} {initial_cluster_profile_consensus} {mmseqs_threads}',

        # mmseqs initial profile search (profile vs consensus), cluster and merger
        'mmseqs_initial_search': f'{exe_mmseqs} search {initial_cluster_profile} {initial_cluster_profile_consensus} {initial_search} {path_mmseqs_tmp} {mmseqs_threads} {extra_opts_search}',
        'mmseqs_clust': f'{exe_mmseqs} clust {initial_cluster_profile_consensus} {initial_search} {search_clust} {mmseqs_threads} {extra_opts_clust}',
        'mmseqs_merge': f'{exe_mmseqs} mergeclusters {seqDB} {search_clust_merged} {initial_cluster} {search_clust}',

        'mmseqs_lndb': f'{exe_mmseqs} lndb {search_clust_merged} {root}/{basename}.profile_cluster',
        'mmseqs_createtsv': f'{exe_mmseqs} createtsv {seqDB} {seqDB} {root}/{basename}.profile_cluster {root}/{basename}.profile_cluster.tsv'

    }
    # execute commands
    for key, value in run_dict.items():
        print(f'{key}: {value}\n')
        subprocess.run(value.split())

    return f'{root}/{basename}.profile_cluster'


# iterate profile and consensus creation, search nad clsutering as per cascaded profile clustering
def cascade_profile_cluster(root, seqDB, clustDB, cascade_steps, threads, extra_opts=None):
    from paths_and_parameters import path_mmseqs_tmp, exe_mmseqs

    # define folder structure
    subprocess.run(f'mkdir {root}/cascaded/'.split())

    # define file names
    basename = seqDB.split('/')[-1]
    casc_clust_merged = f'{root}/cascaded/{basename}.casc-clust'
    all_casc_cluster_steps = []

    # define and assign commands and parameters to run
    mmseqs_threads = f'--threads {threads}'
    extra_opts_cascade_cluster = ''
    extra_opts_profile = ''
    extra_opts_search = ''
    extra_opts_clust = ''

    if isinstance(extra_opts, dict):
        try:
            extra_opts_profile = extra_opts['result2profile']
            extra_opts_search = extra_opts['cascade_search']
            extra_opts_clust = extra_opts['cascade_clust']

        except KeyError:
            print(f'WARNING: Extra opts specified for initial_profile_cluster() but not all found in keys!')

    # iterate over number of steps, create a folder for each step
    # starting with cluster result from initial_profile_cluster()

    iter_clustDB = clustDB

    for iter in range(1, cascade_steps + 1):

        iter_out = f'{root}/cascaded/cascade_{iter}/'
        iter_basename = f'{iter_out}{basename}.casc_{iter}'
        casc_profiles = f'{iter_basename}.profiles'
        casc_consensus = f'{iter_basename}.consensus'
        casc_search = f'{iter_basename}.search'
        casc_clust = f'{iter_basename}.search.clust'

        # read
        run_dict = {'mkdir': f'mkdir {iter_out}',
                    'mmseqs_casc_r2p': f'{exe_mmseqs} result2profile {seqDB} {seqDB} {iter_clustDB} {casc_profiles} {mmseqs_threads} {extra_opts_profile}',
                    'mmseqs_casc_p2c': f'{exe_mmseqs} profile2consensus {casc_profiles} {casc_consensus} {mmseqs_threads}',
                    'mmseqs_casc_search': f'{exe_mmseqs} search {casc_profiles} {casc_consensus} {casc_search} {path_mmseqs_tmp} {mmseqs_threads} {extra_opts_search}',
                    'mmseqs_casc_clust': f'{exe_mmseqs} clust {casc_consensus} {casc_search} {casc_clust} {mmseqs_threads} {extra_opts_clust}'

                    }
        # execute commands
        for key, value in run_dict.items():
            print(f'{key}: {value}\n')
            subprocess.run(value.split())

        # update clustDB for next loop
        iter_clustDB = casc_clust

        # store final output clusters for last merger
        all_casc_cluster_steps.append(casc_clust)

    # merge cascaded clusters and link final result
    mmseqs_casc_merge = f'{exe_mmseqs} mergeclusters {seqDB} {casc_clust_merged} {clustDB} {" ".join(all_casc_cluster_steps)}'
    mmseqs_lndb = f'{exe_mmseqs} lndb {casc_clust_merged} {root}/{basename}.cascaded_cluster'
    mmseqs_createtsv = f'{exe_mmseqs} createtsv {seqDB} {seqDB} {root}/{basename}.cascaded_cluster {root}/{basename}.cascaded_cluster.tsv'
    subprocess.run(mmseqs_casc_merge.split())
    subprocess.run(mmseqs_lndb.split())
    subprocess.run(mmseqs_createtsv.split())

    return f'{root}/{basename}.cascaded_cluster'


# write fasta files into individual files for realignment
def extract_query_fasta(root, seqDB, clusterDB, skip_singletons=True):
    from paths_and_parameters import exe_mmseqs

    # define and assign commands and parameters to run
    basename = clusterDB.split('/')[-1]
    fasta_root = f'{root}cluster_fastas/'
    fastaDB = f'{root}{basename}.fasta'

    print(basename)
    print(fasta_root)
    print(fastaDB)

    # createseqdb to get cluster members as \x00 delimited fasta
    mkdir = f'mkdir {fasta_root}'
    mmseqs_subdb = f'{exe_mmseqs} createseqfiledb {seqDB} {clusterDB} {fastaDB}'
    subprocess.run(mkdir.split())
    subprocess.run(mmseqs_subdb.split())

    # read the binary mmseqsDB and split at \x00 bytes. Write each part into individual file.
    with open(fastaDB, 'rb') as binary_fasta:
        # drop last empty slice
        cluster_fastas = binary_fasta.read().split(b'\x00')[:-1]

    for fasta in cluster_fastas:

        fastastr = fasta.decode('utf-8', 'ignore')
        name = fastastr.split('\n')[0].strip('>').split(' ')[0]
        size = fastastr.count('>')

        if skip_singletons and size > 1:
            with open(f'{fasta_root}{name}', 'w') as fastafile:
                fastafile.write(fastastr)

        elif not skip_singletons:
            with open(f'{fasta_root}{name}', 'w') as fastafile:
                fastafile.write(fastastr)

        # print(f'Wrote {n} {fasta_root}{name}')

    return fasta_root


# direct realign all fasta files in directory
def realign_all_fastas(fasta_root, threads, max_sequences=20, save_intermediate_files=False):
    from muscle_crop_and_align import realign_cluster_fasta

    # get all fasta files in folder
    ls_command = f'ls {fasta_root}*.fasta'
    fastafiles = subprocess.run(ls_command.split(), capture_output=True, text=True).stdout.split()

    for fasta in fastafiles:
        fastafile = fasta_root + fasta
        realign_cluster_fasta(fastafile, threads, max_sequences=20, save_intermediate_files=False)

    return


# realign all fasta files in directory submitting individual batch jobs to the swarm system
# swarm header formatted from swarm opts dictionary of parameters
# requires ONLY fasta files in fasta_root
# this is as we avoid .fasta extension as hhmake takes index name from file name
#OBSOLETE!!!
def realign_all_fastas_swarm(fasta_root, swarm_opts, max_sequences=250, filter_entropy=0.5, muscle_reps=10):
    from paths_and_parameters import path_tmp, exe_python, exe_realignment
    from core_functions.helper_functions import swarm_submit_with_lock

    swarmfile = path_tmp + 'realign_all_fastas.swarm'

    # get all files in folder
    ls_command = f"find {fasta_root} -type f"
    fastafiles = subprocess.run(ls_command.split(), capture_output=True, text=True).stdout
    fastafiles = [file for file in fastafiles.split('\n') if file != '']

    swarm = open(swarmfile, 'w')

    # format swarm header for submission
    threads = swarm_opts['threads-per-process']
    swarm.write(''.join([f'#SWARM --{key} {value}\n' for key, value in swarm_opts.items()]))

    for fasta in fastafiles:

        # exclude all singletons from realignment submission
        with open(fasta, 'r') as fastafile:
            size = fastafile.read().count('>')

        if size > 1:
            submit_command_string = f'{exe_python} -u {exe_realignment} --fasta {fasta} --threads {threads} --max_seqs {max_sequences} --filter_entropy {filter_entropy} --muscle_reps {muscle_reps}'
            swarm.write(submit_command_string + '\n')

    swarm.close()

    swarm_submit_with_lock(swarmfile, refresh=60)

    return fasta_root


# realign all fasta files in directory submitting individual batch jobs to the swarm system
# swarm header formatted from swarm opts dictionary of parameters
# requires ONLY fasta files in fasta_root
# this is as we avoid .fasta extension as hhmake takes index name from file name
def realign_all_fastas_swarm_binned(fasta_root, swarm_opts, max_swarms=1000, max_sequences=250, filter_entropy=0.5, muscle_reps=10):
    from paths_and_parameters import path_tmp, exe_python, exe_realignment, muscle_realignment_timeout
    from core_functions.helper_functions import swarm_submit_with_lock, greedy_pack

    swarmfile = path_tmp + 'realign_all_fastas.swarm'

    # get all files in folder
    ls_command = f"find {fasta_root} -type f"
    fastafiles = subprocess.run(ls_command.split(), capture_output=True, text=True).stdout
    fastafiles = [file for file in fastafiles.split('\n') if file != '']

    swarm = open(swarmfile, 'w')

    # format swarm header for submission
    threads = swarm_opts['threads-per-process']
    swarm.write(''.join([f'#SWARM --{key} {value}\n' for key, value in swarm_opts.items()]))

    #format the item dictionary for bin packing
    items_to_pack = []

    for fasta in fastafiles:
        #calculate the size of the file for batching
        with open(fasta, 'r') as fastafile:
            fastastr = fastafile.read()
            size = fastastr.count('>')

            # dont realign singletons
            if size > 1:
                # ad_hoc estimate for the runnning time complexity as O((mean_len*num)**1.4) adjusted for max_alignment size
                weight = int(((len(fastastr) / size) * min([size, max_sequences])) ** 1.4)
                items_to_pack.append((fasta, weight))


    packed_items, packed_weights = greedy_pack(items_to_pack, max_swarms)
    print(f'Packed {len(items_to_pack)} files into {len(packed_items)} bins')

    for bin in packed_items:
        command_list = []
        for fasta in bin:
            command_list.append(f'{exe_python} -u {exe_realignment} --fasta {fasta} --threads {threads} --max_seqs {max_sequences} --filter_entropy {filter_entropy} --muscle_reps {muscle_reps} --muscle_timeout {muscle_realignment_timeout}')

        swarm.write('; '.join(command_list) + '\n')

    swarm.close()

    swarm_submit_with_lock(swarmfile, refresh=60)

    return fasta_root


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
