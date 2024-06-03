import argparse

from core_functions.structure_alignment import run_ICARUS

parser = argparse.ArgumentParser(description='Reduce fasta to size and realign using multiple rounds of muscle')
parser.add_argument('--root', type=str, required=True, help='output root')
parser.add_argument('--pdb1', type=str, required=True, help='moving pdb absolute path')
parser.add_argument('--pdb2', type=str, required=True, help='reference pdb absolute path')
parser.add_argument('--level', type=int, required=True, help='ICARUS domain partitioning level')
parser.add_argument('--threads', type=int, required=True, help='threads to run for all subprocesses')
parser.add_argument('--save_intermediate_files', action=argparse.BooleanOptionalAction)
args = parser.parse_args()

#run main
if __name__ == '__main__':

    result_name = run_ICARUS(root = args.root,
                             pdb1 = args.pdb1,
                             pdb2 = args.pdb2,
                             level = args.level,
                             threads = args.threads,
                             #save_intermediate_files=args.save_intermediate_files
                             )