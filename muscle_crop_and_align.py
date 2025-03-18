import argparse
import subprocess


from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy
from core_functions.microcosm_functions import fasta_reduce_size
from core_functions.software_wrappers import muscle_ensamble

#realign target fastafile for profile generation with fasttree leaf-cropping for large fastas
def realign_cluster_fasta(fasta, threads, max_sequences=250, filter_entropy=0, gap_fraction=0.5, muscle_reps=10, save_intermediate_files=False, muscle_timeout=7200):

    #configure paths
    basename = fasta.split('/')[-1]

    with open(fasta, 'r') as fastafile:
        size = fastafile.read().count('>')

    if size == 1:
        print(f' There is only one seuqnce in {fasta} ({size}), will skip')
        return

    # if there are too many sequences, crop to size
    if size > max_sequences:
        print(f' There are more than {max_sequences} sequences in {fasta} ({size}), will crop to size')

        #filter entropy for cropping always 0
        crop_entropy = 0

        fasta = fasta_reduce_size(fasta, threads, max_sequences, crop_entropy, save_intermediate_files=save_intermediate_files)

    else:
        print(f' There are less than {max_sequences} sequences in {fasta} ({size}), no cropping needed')

    # align
    print(f' Aligning with muscle5 as ensemble with {muscle_reps} replicates')
    muscle_fasta = muscle_ensamble(fasta, threads, muscle_reps, super5=False, save_intermediate_files=save_intermediate_files, timeout=muscle_timeout)

    # filter by entropy
    aln = fasta_to_dict(file=muscle_fasta)
    aln_filter = filter_by_entropy(aln, filter_entropy, gap_fraction)

    if save_intermediate_files:
        dict_to_fasta(aln_filter, write_file=f'{fasta}.aln')

    else:
        dict_to_fasta(aln_filter, write_file=fasta)
        subprocess.run(f'rm {fasta}.muscle {fasta}.muscle-efa {fasta}.muscle.log'.split())

    return

#argparse define
parser = argparse.ArgumentParser(description='Reduce fasta to size and realign using multiple rounds of muscle')
parser.add_argument('--fasta', type=str, required=True, help='fasta file')
parser.add_argument('--threads', type=int, required=False, nargs='?', const=1, default=1, help='threads to run for all subprocesses')
parser.add_argument('--max_seqs', type=int, required=False, nargs='?', const=1, default=300, help='maximum allowed sequences without filtering')
parser.add_argument('--filter_entropy', type=float, required=False, nargs='?', const=1, default=0, help='minimal column entropy to allow in alignment')
parser.add_argument('--gap_fraction', type=float, required=False, nargs='?', const=1, default=1, help='fraction of total gaps which to allow, ex 0.2 allows maximum 20 gaps per column')
parser.add_argument('--muscle_reps', type=int, required=False, nargs='?', const=1, default=3, help='number of diversified muscle iterations, must be greater than 1')
parser.add_argument('--muscle_timeout', type=float, required=False, nargs='?', const=1, default=3600, help='soft running time timeout, if exceeded will resubmit with super5')
parser.add_argument('--save_intermediate_files', action=argparse.BooleanOptionalAction)
args = parser.parse_args()

#run main
if __name__ == '__main__':

    realign_cluster_fasta(fasta=args.fasta,
                          threads=args.threads,
                          max_sequences=args.max_seqs,
                          filter_entropy=args.filter_entropy,
                          gap_fraction=args.gap_fraction,
                          muscle_reps=args.muscle_reps,
                          save_intermediate_files=args.save_intermediate_files,
                          muscle_timeout=args.muscle_timeout)
