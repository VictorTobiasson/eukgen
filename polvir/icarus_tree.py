import argparse
import sys
import os
import subprocess
import numpy as np
module_path = "/data/luojaa/eukgen/"
sys.path.append(module_path)

#prepare all mmseqs
def icarus_tree_run(icarus_root, alignment_root, alignment_basename, threads, save_intermediate_files):
    
    from paths_and_parameters import microcosm_format_opts
    from core_functions.software_wrappers import calculate_IQtree
    from polvir.parseIcarus import parse_alignment_summary, define_gaps, map_gaps, crop_pdb, entropy_filter
    # # run af2 script
    # # run icarus
    
    # # parse icarus output and crop PDB files
    # for seq in os.listdir(icarus_root):
    #     pdb_renum = f"{icarus_root}{seq}/icarus_results_on_ED_38727.Pavir.Eb01025.1/models/{seq}.renum.pdb"
    #     struct_align = f"{icarus_root}{seq}/icarus_results_on_ED_38727.Pavir.Eb01025.1/.icarus*/icarus_output/results/PDB1_and_PDB2/summary.txt"
    #     alignment_data = parse_alignment_summary(struct_align)
    #     gaps = define_gaps(alignment_data)
    #     gap_ranges = map_gaps(gaps)
    #     crop_pdb(icarus_root, alignment_root, seq, gap_ranges)

    # # create output fasta files and run alignments in sbatch
    # dists = alignment_data["dist"]
    # query = alignment_data["QUERY"]
    # dists_qinserts = "".join(list(np.array(list(dists))[(np.array(list(query)) != "-") & (np.array(list(query)) != " ")]))
    # query_contig = query.replace(" ", "").replace("-","")
    # notgaps = [not gap for gap in gaps]
    # dists_clean = "".join(list(np.array(list(dists_qinserts))[notgaps]))
    # query_clean = "".join(list(np.array(list(query_contig))[notgaps]))
    # outfasta = f"{alignment_root}input_fastas/{seq}.fasta"
    # with open(outfasta,"w") as outfile:
    #     print(f">{seq}", file=outfile)
    #     print(query_clean, file=outfile)
    # concatfasta = f"{alignment_root}icarus_cropped.fasta"
    # subprocess.run(f"cat {alignment_root}input_fastas/*.fasta > {concatfasta}", shell=True)
    
    # entropy filter
    entropy_filter(alignment_root)
    
    # iqtree write swarm file

    # run IQtree
    align_fasta = f"{alignment_root}{alignment_basename}.entropy_filt.fasta"
    calculate_IQtree(align_fasta,
                    evo_model_params = microcosm_format_opts['evo_model_params'],
                    threads = threads,
                    save_intermediate_files=save_intermediate_files)


#argparse define
parser = argparse.ArgumentParser(description='Run IQTree on Icarus informed alignments')
parser.add_argument('--icarus_root', type=str, required=True, help='path to icarus output root')
parser.add_argument('--alignment_root', type=str, required=True, help='path to icarus alignment root')
parser.add_argument('--alignment_basename', type=str, required=True, help='alignment file basename, ex: "poly_380.foldmason_aa", "icarus_280.mafft_linsi"')
parser.add_argument('--threads', type=int, nargs='?', default=-1, const=-1, help='number of threads')
parser.add_argument('--delete_intermediate_files', action='store_false', help='remove all intermediate files')

args = parser.parse_args()

#run main
if __name__ == '__main__':

    icarus_tree_run(args.icarus_root,
                    args.alignment_root,
                    args.alignment_basename,
                    args.threads,
                    save_intermediate_files=args.delete_intermediate_files)
