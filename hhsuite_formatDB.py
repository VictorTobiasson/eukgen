import argparse
import subprocess


from core_functions.hhsuite_functions import build_HHsuite_database_from_MSAs



#argparse define
parser = argparse.ArgumentParser(description='Reduce fasta to size and realign using multiple rounds of muscle')
parser.add_argument('--fasta_root', type=str, required=True, help='folder of MSAs, should ONLY contain MSAS')
parser.add_argument('--output_root', type=str, required=True, help='folder for output DB')
parser.add_argument('--db_name', type=str, required=True, help='root_name for _hhm, _cs219 and _a3m files')
parser.add_argument('--mpis', type=int, required=True, help='parallel mpis')

args = parser.parse_args()

#run main
if __name__ == '__main__':

    build_HHsuite_database_from_MSAs(fasta_root = args.fasta_root,
                                     output_root = args.output_root,
                                     DB_name = args.db_name,
                                     mpis = args.mpis)
