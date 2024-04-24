from multiprocessing import current_process
import subprocess


# run command for generating a diversified ensemble with muscle and then extracting the maxcc aln
def muscle_ensamble(base_fasta, threads, muscle_reps, super5=False, save_intermediate_files=True, timeout=None):

    thread = current_process().pid
    threadID_string = f'{thread} | {base_fasta}:'

    from paths_and_parameters import exe_muscle5
    import os

    os.environ['OMP_NUM_THREADS'] = str(threads)

    logfile = open(f'{base_fasta}.muscle.log', 'a')

    if super5:
        for rep in range(muscle_reps):
            print(f'Aligning super5 iteration {rep}')
            muscle5_align_ete_params = f" -super5 {base_fasta} -output {base_fasta}.@.super5-tmp -perm all -perturb {rep} -threads {threads}"
            muscle5_ete_command = exe_muscle5 + muscle5_align_ete_params

            # run with per iteration timeout, rerun with single muscle_rep if timeout still exceeded
            try:
                subprocess.run(muscle5_ete_command.split(), stdout=logfile, stderr=logfile, timeout=timeout)

            except subprocess.TimeoutExpired:
                print(f'Iteration timeout exceeded for {base_fasta} rerunning super5 with 1 replicate')

                #rerun one perturbed super5 replicate
                muscle_ensamble(base_fasta, threads, muscle_reps=1, super5=True, save_intermediate_files=save_intermediate_files, timeout=None)

                return f'{base_fasta}.muscle'

        #extract folder name
        full_file_path = subprocess.run(f'realpath {base_fasta}'.split(), capture_output=True, text=True).stdout.strip()
        root = '/'.join(full_file_path.split('/')[:-1]) + '/'
        name = base_fasta.split('/')[-1]

        ls_command = f"find {root} -name '{name}*super5-tmp' -type f"
        efa_files = subprocess.run(ls_command, shell=True, capture_output=True, text=True).stdout
        efa_files = [file for file in efa_files.split('\n') if file != '']

        with open(base_fasta + '.muscle-efa', 'w') as efa_merge:
            for file in efa_files:
                with open(file, 'r') as efa_in:
                    efa_merge.write(f'<{file}\n')
                    efa_merge.write(efa_in.read())

                    if not save_intermediate_files:
                        subprocess.run(f'rm {file}'.split())


    else:
        # run diversified ensemble
        muscle5_align_ete_params = f" -threads {threads} -diversified -replicates {muscle_reps} -align {base_fasta} -output {base_fasta}.muscle-efa"
        muscle5_ete_command = exe_muscle5 + muscle5_align_ete_params

        #run with global timeout, submit as super5 if timeout exceeded
        try:
            subprocess.run(muscle5_ete_command.split(), stdout=logfile, stderr=logfile, timeout=timeout)
            print(muscle5_ete_command)

            print(muscle_reps, super5)
            if muscle_reps == 1 and not super5:
                # mv .muscle-efa to .muscle for single iteration runs
                subprocess.run(f'mv {base_fasta}.muscle-efa {base_fasta}.muscle'.split(), stdout=logfile, stderr=logfile)

        except subprocess.TimeoutExpired:
            print(f'Timeout exceeded for {base_fasta} rerunning with super5')
            muscle_ensamble(base_fasta, threads, muscle_reps, super5=True, save_intermediate_files=save_intermediate_files, timeout=timeout)

    if super5 or muscle_reps > 1:
        # for multiple repetitions select -maxcc
        print('maxxcc merging')
        muscle5_maxcc_params = f' -maxcc {base_fasta}.muscle-efa -output {base_fasta}.muscle'
        muscle5_command = exe_muscle5 + muscle5_maxcc_params
        subprocess.run(muscle5_command.split(), stdout=logfile, stderr=logfile)

    logfile.close()

    return f'{base_fasta}.muscle'

# basic wrapper to launch IQTree2 and save output and errors to separate file
def calculate_IQtree(alignment, evo_model_params, threads=1, save_intermediate_files=True):

    from paths_and_parameters import exe_iqtree

    # configure paths and flags
    thread = current_process().pid
    threadID_string = f'{thread} | :'

    # construct IQtree
    print(threadID_string + f' Constructing IQTree for {alignment}')

    with open(f'{alignment}.iqtree.log', 'a') as iqtree_logfile:
        # run 1000 ultrafast bootstraps estimating a subset fo models, add -bnni for UFBoot model violations
        iqtree_command = f'{exe_iqtree} -s {alignment} {evo_model_params} --threads {threads} -B 1000 -bnni --redo'
        subprocess.run(iqtree_command.split(), stdout=iqtree_logfile, stderr=iqtree_logfile)

    #save only .treefile
    if not save_intermediate_files:
        subprocess.run(f'rm {alignment}{{.ckp.gz,.contree,.iqtree,.iqtree.log,.splits.nex,.model.gz,.mldist,.bionj}}'.split())

    return f'{alignment}.treefile'