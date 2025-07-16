import argparse
import subprocess


from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy

#realign target fastafile for profile generation with fasttree leaf-cropping for large fastas
def filter_msa(fasta, entropy_min=0, gap_fraction=0.5, gaptoken='-', verbose=True):

    # filter by entropy
    aln = fasta_to_dict(file=fasta)
    
    initial_length = len(aln[list(aln.keys()).pop()])
    
    aln_filter = filter_by_entropy(aln, entropy_min, gap_fraction, seq_length_frac_min=None, filter_accs=[], gaptoken='-')
    dict_to_fasta(aln_filter, write_file=f'{fasta}.filtered', verbose=False)
    
    final_length = len(aln_filter[list(aln_filter.keys()).pop()])
    
    gap_pct = (initial_length-final_length)*100/initial_length
    if verbose:
        print(f'Filtered {fasta} at gap={gap_fraction} and entropy={entropy_min}')
        print(f'Filter removed {initial_length-final_length} columns of {initial_length} ({round(gap_pct, 1)}%)')
        print(f'Wrote output {fasta}.filtered')


# +
#argparse define
parser = argparse.ArgumentParser(description='Filter fasta msa file by gaps and/or entropy')
parser.add_argument('--fasta', type=str, required=True, help='Input fasta file')
parser.add_argument('--entropy_min', type=float, required=False, nargs='?', const=1, default=0, help='minimal column entropy to allow in alignment')
parser.add_argument('--gap_fraction', type=float, required=False, nargs='?', const=1, default=1, help='fraction of total gaps which to allow, ex 0.2 allows maximum 20%% gaps per column')
parser.add_argument('--gap_token', type=str, required=False, nargs='?', const='-', help='token denominating gaps in fasta MSA')
#parser.add_argument('--verbose', type=bool, required=False, nargs='?', const=False, help='output gap statistics')
parser.add_argument('--verbose', action=argparse.BooleanOptionalAction)

args = parser.parse_args()
# -

#run main
if __name__ == '__main__':
    
    filter_msa(fasta=args.fasta, 
               entropy_min=args.entropy_min, 
               gap_fraction=args.gap_fraction,
               gaptoken=args.gap_token,
               verbose=args.verbose)
