import os
import argparse
import subprocess
import sys
module_path = "/data/luojaa/eukgen/"
if module_path not in sys.path:
    sys.path.append(module_path)

def merge_hhsuite_search_results_chunks(search_root, output_basename, chunks, pairwise_cov = False, probability = False):
    from core_functions.hhsuite_functions import merge_hhsuite_search_results
    from core_functions.helper_functions import swarm_submit_with_lock
    from paths_and_parameters import exe_python, exe_hhsuite_parse, parse_swarm_opts
    
    #create chunks first
    print(f"creating {chunks} chunks")
    results_files = os.listdir(search_root)
    results_paths = [search_root + f for f in results_files]
    numfiles = len(results_paths)
    chunksize = numfiles // chunks
    for chunk_no in range(chunks): # CHANGE TO ALL CHUNKS AFTER TESTING
        start = chunksize * chunk_no
        if chunk_no == chunks - 1:
            end = numfiles
        else:
            end = chunksize * (chunk_no + 1)
        chunk_files = results_paths[start:end]
        chunkdir = f"{search_root}chunk{chunk_no}/"
        process = subprocess.Popen(["mkdir", chunkdir[:-1]])
        process.wait()
        for file in chunk_files:
            process = subprocess.Popen(["mv", file, chunkdir]) # CHANGE TO MV
            process.wait()
    
    #probably better to swarm this
    # format swarm file
    print("formatting swarm file")
    os.system(f"mkdir {search_root}tmp/")
    swarmfile = search_root + 'tmp/hhsuite_format.swarm'
    swarm = open(swarmfile, 'w')

    # format swarm header for submission
    swarm.write(''.join([f'#SWARM --{key} {value}\n' for key, value in parse_swarm_opts.items()]))
    
    for chunk_no in range(chunks):
        #add execution line to swarmfile
        swarm.write(f'{exe_python} -u {exe_hhsuite_parse} --search_root {search_root} --output_basename {output_basename} --chunk_no {chunk_no} --pairwise_cov {pairwise_cov} --probability {probability} \n')

    swarm.close()

    swarm_submit_with_lock(swarmfile, refresh=60)

    


#argparse define
parser = argparse.ArgumentParser(description='parse and filter hhblits search results, serially in chunks to reduce memory use')
parser.add_argument('--search_root', type=str, required=True, help='input path to hhblits results </data/user/.../hhsuite/>XXX_hmm.{ffindex,ffdata}')
parser.add_argument('--output_basename', type=str, required=True, help='output path for merged hhblits output </data/user/.../YYY>.{tsv,query.tsv,pkl}')
parser.add_argument('--chunks', type=int, required=True, help='number of chunks')
parser.add_argument('--pairwise_cov', type=float, required=False, nargs='?', const=False, help='pairwise coverage filter for hhblits results [0,1]')
parser.add_argument('--probability', type=int, required=False, nargs='?', const=False, help='probability filter for hhblits results [0, 100]')
args = parser.parse_args()

#run main
if __name__ == '__main__':
    merge_hhsuite_search_results_chunks(search_root=args.search_root,
                   output_basename = args.output_basename,
                   chunks = args.chunks,
                   pairwise_cov = args.pairwise_cov,
                   probability = args.probability,
                   )
