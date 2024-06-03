# python global packages
import pandas as pd
import subprocess
import os


def get_foldseek_mediod(root, basename, pdbdir, threads=16, delete_intermediate_files=False):

    from paths_and_parameters import exe_foldseek

    search_opts = f"-s 2 --format-mode 4 --format-output 'query,target,alntmscore' --threads {threads}"

    run_dict = {
        'foldseek_tmp_rundir': f'mkdir {root}/foldseek_mediod/',
        'foldseek_easy-search': f'{exe_foldseek} easy-search {pdbdir} {pdbdir} {root}/foldseek_mediod/foldseek_tmscores {root}/foldseek_mediod/ {search_opts}',
        'foldseek_cp_results': f'cp {root}/foldseek_mediod/foldseek_tmscores {root}/{basename}.tmscores.tsv'
    }

    # execute commands
    for key, value in run_dict.items():
        print(f'{key}: {value}\n')
        subprocess.run(value, shell=True)

    # optional delete
    if delete_intermediate_files:
        subprocess.run(f'rm -r {root}/foldseek_mediod/', shell=True)

    # get mediod
    tmdata = pd.read_csv(f'{root}/{basename}.tmscores.tsv', sep='\t', index_col=0)
    mediod = tmdata.groupby('query').apply(lambda x: x['alntmscore'].sum()).idxmin()

    return mediod, f'{root}/{basename}.tmscores.tsv'


# very messy requiring a temporary directory shift to execute icarus as it always writes inot current directory.
def run_ICARUS(root, pdb1, pdb2, threads, level=1):

    from paths_and_parameters import exe_icarus
    from uuid import uuid4

    # icarus path config
    pdb1_name = pdb1.split('/')[-1]
    pdb1_base = pdb1_name.strip(".pdb")
    pdb2_name = pdb2.split('/')[-1]
    pdb2_base = pdb2_name.strip(".pdb")

    align_opts = f'-l {level} -c {threads} --verbose --force'
    icarus_output_root = f'icarus_output/results/PDB1_and_PDB2/PDB1_on_PDB2/'
    icarus_output_models = f'{icarus_output_root}/result_PDBs/'
    icarus_output_data = f'{icarus_output_root}/gdt*'

    # icarus is sensitive to running directory, we create a tmp and final merged output directory for each run
    # icarus also does not handle filenames well so rename the temp files to PDB1
    output_root = f'icarus_results_on_{pdb2_base}'
    origin_dir = subprocess.run('pwd', capture_output=True, text=True).stdout.strip('\n')
    tmp_dir = f'.icarus_tmp.{uuid4().hex}'

    os.makedirs(f'{root}/', exist_ok=True)
    os.makedirs(f'{root}/{output_root}/', exist_ok=True)
    os.makedirs(f'{root}/{output_root}/models', exist_ok=True)
    os.makedirs(f'{root}/{output_root}/data', exist_ok=True)
    os.makedirs(f'{root}/{output_root}/{tmp_dir}', exist_ok=True)
    os.chdir(f'{root}/{output_root}/{tmp_dir}')

    print(f'{root}/')
    print(f'{root}/{output_root}/')
    print(f'{root}/{output_root}/{tmp_dir}')

    # define commands, run in tmp_dir
    run_dict = {
        'icarus_cp_input1': f'cp {pdb1} ./PDB1.pdb',
        'icarus_cp_input2': f'cp {pdb2} ./PDB2.pdb',
        'icarus_align': f'{exe_icarus} -p1 PDB1.pdb -p2 PDB2.pdb {align_opts}',
        'icarus_mv_data': f'mv {icarus_output_data} ../data/{pdb1_base}.icarus'

    }

    # execute commands
    for key, value in run_dict.items():
        print(f'\n{key}: {value}\n')
        subprocess.run(value, shell=True)

    print('DONE')

    # icarus does not store or provide information for what is the best model relative to one reference
    # best model is most often the latest written.
    result_files = subprocess.run(f'find {icarus_output_models} -type f'.split(), capture_output=True, text=True).stdout
    pdb_files = [file for file in result_files.split('\n') if file.split('PUs')[-1] in ['.pdb']]
    latest_file = sorted(pdb_files, key=os.path.getctime)[0]
    subprocess.run(f'cp {latest_file} ../models/{pdb1_name}'.split())

    # return to original directory and remove tmp
    os.chdir(origin_dir)
    subprocess.run(f'rm -rf {tmp_dir}'.split())

    return output_root

def format_ICARUS_swarm(root, pdb_root, target, threads, level, swarm_opts):

    from paths_and_parameters import path_tmp, exe_python, exe_icarus_alignment, exe_icarus_environment
    from core_functions.helper_functions import swarm_submit_with_lock

    swarmfile = path_tmp + 'icarus_alignment.swarm'

    # get all files in folder
    ls_command = f"find {pdb_root} -type f"
    pdbfiles = subprocess.run(ls_command.split(), capture_output=True, text=True).stdout
    pdbfiles = [file for file in pdbfiles.split('\n') if file not in ['', target]]

    # format swarm header for submission
    swarm = open(swarmfile, 'w')
    swarm.write(''.join([f'#SWARM --{key} {value}\n' for key, value in swarm_opts.items()]))

    print(f'Will align {len(pdbfiles)} structures against {target}')

    for pdb in pdbfiles:
        submit_conda_string = f'source myconda; conda activate {exe_icarus_environment}'
        submit_command_string = f'{exe_python} -u {exe_icarus_alignment} --root {root} --pdb1 {pdb} --pdb2 {target} --level {level} --threads {threads} --save_intermediate_files'
        swarm.write(submit_conda_string + '; ' + submit_command_string + '\n')

    swarm.close()

    swarm_submit_with_lock(swarmfile, refresh=60)

    return

