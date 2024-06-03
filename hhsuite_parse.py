import os
import argparse
import subprocess
import sys

module_path = "/data/luojaa/eukgen/"
if module_path not in sys.path:
    sys.path.append(module_path)

def merge_hhsuite_search_results_exe(search_root, output_basename, chunk_no, pairwise_cov = False, probability = False):
    from core_functions.hhsuite_functions import merge_hhsuite_search_results

    chunkdir = f"{search_root}chunk{chunk_no}/"
    output_name = f"{output_basename}chunk{chunk_no}"
    merge_hhsuite_search_results(chunkdir, output_name, write_tsv=True, filter_cov=pairwise_cov, filter_prob=probability)


#argparse define
parser = argparse.ArgumentParser(description='parse and filter hhblits search results, serially in chunks to reduce memory use')
parser.add_argument('--search_root', type=str, required=True, help='input path to hhblits results </data/user/.../hhsuite/>XXX_hmm.{ffindex,ffdata}')
parser.add_argument('--output_basename', type=str, required=True, help='output path for merged hhblits output </data/user/.../YYY>.{tsv,query.tsv,pkl}')
parser.add_argument('--chunk_no', type=int, required=True, help='chunk to parse')
parser.add_argument('--pairwise_cov', type=float, required=False, nargs='?', const=False, help='pairwise coverage filter for hhblits results [0,1]')
parser.add_argument('--probability', type=int, required=False, nargs='?', const=False, help='probability filter for hhblits results [0, 100]')
args = parser.parse_args()

#run main
if __name__ == '__main__':
    merge_hhsuite_search_results_exe(search_root=args.search_root,
                   output_basename = args.output_basename,
                   chunk_no = args.chunk_no,
                   pairwise_cov = args.pairwise_cov,
                   probability = args.probability,
                   )
