import ete3
from Bio import SeqIO
import regex as re
import os

base = "/data/luojaa/kegg/kog_uid/"
kog2uidpaths = [base + tsv for tsv in os.listdir(base)]
uid2kog = {} 
for path in kog2uidpaths:
    with open(path, "r") as f:
        for line in f:
            kog, uid, descriptor = line.strip().split("\t")
            uid2kog[uid] = kog
print("uid2kogs mapped")            

def extract_taxa(instr):
    taxa = re.findall(r"OX=([0-9]+)\s", instr)[0]
    return taxa
def extract_uid(instr):
    uid = instr.split("|")[1]
    return uid

ncbi = ete3.NCBITaxa()
def isEuk(taxid):
    lineage_ids = ncbi.get_lineage(taxid)
    if 2759 in lineage_ids:
        return True
    return False

in_fasta = "/data/luojaa/uniprot_combined.fasta"
uid_taxids = "/data/luojaa/uid2taxids2iseuk.combined.tsv"
log_path = "/data/luojaa/log/parse_fastas.combined.prok.log"
skip = True

records_count, log_count, euk_count = 0, 0, 0
with open(in_fasta, "r") as handle, open(uid_taxids, "a") as outfile:
    # iterate through trembl fasta file: 200M proteins
    for record in SeqIO.parse(handle, 'fasta'):
        if records_count % 100000 == 0:
            with open(log_path, "a") as logfile:
                print(",".join([str(records_count), "records_considered"]), file=logfile)
        records_count += 1
        identifier, description, sequence = record.id, record.description, record.seq
        uid = extract_uid(identifier)
        taxid = extract_taxa(description)
        # # skip all entries until fail point
        # if (taxid == '1554493') or (taxid == 1554493):
        #     skip = False
        # if skip:
        #     continue
        # does this trembl entry appear in our list of kog protein uids?
        try:
            kog = uid2kog[uid]
            if log_count % 10000 == 0:
                with open(log_path, "a") as logfile:
                    print(",".join([str(log_count), "uids_mapped"]), file=logfile)
            log_count +=1
        # if not, go next
        except KeyError:
            continue
        # if so, is this entry a eukaryotic gene?
        try:
            if not isEuk(taxid):
                print("\t".join([uid, taxid, "False"]), file=outfile)
                if euk_count % 1000 == 0:
                    with open(log_path, "a") as logfile:
                        print(",".join([str(euk_count), "proks_recorded"]), file=logfile)
                # add this gene to it's kog's fasta file
                outfasta=f"/data/luojaa/kog_fastas_prok/{kog}.fasta"
                with open(outfasta, "a") as cluster:
                    print(">" + uid, file=cluster)
                    print(sequence, file=cluster)
                euk_count += 1
            else:
                print("\t".join([uid, taxid, "True"]), file=outfile)
        except:
            ete3logfile = "/data/luojaa/log/ete3log.prok.fasta"
            with open(ete3logfile, "a") as logfile:
                print(" ".join([kog, uid, taxid, description, "taxid not found"]), file=logfile)
            continue
        
