import subprocess
import regex as re
import ast
import os
from collections import defaultdict

# map kegg ids to ncbi ids
def read_ids():
    idpairs = defaultdict(list)
    base_dir = "Data/Kegg2Uniprot/"
    id_tsvs = [base_dir + tsv for tsv in os.listdir(base_dir) if ("chunk" in tsv) and (".gz" not in tsv)]
    for tsv in id_tsvs:
        with open(tsv, "r") as f:
            next(f)
            lines = f.readlines()
            for line in lines:
                keggid, uniprotid = line.strip().split("\t")
                idpairs[uniprotid].append([keggid])
    return idpairs

# uniprot_ids = "\t".join(list(read_ids().keys())) # saved to Data/uniprot_ids.txt
# with open("Data/uniprot_ids.tsv", "a") as outfile:
#     print(uniprot_ids, file=outfile)


## run batch requests on UNIPROT ids

def run_batch(ids):
    requesturl = "https://rest.uniprot.org/idmapping/run"
    from_db = "\'from=\"UniProtKB_AC-ID\"\'"
    to_db = "\'to=\"RefSeq_Protein\"\'"
    idstr = "\'ids=\"" + ",".join(ids) + "\"\'"
    command = "curl --form {} --form {} --form {} {}".format(from_db, to_db, idstr, requesturl)
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
        shell = True)
    dict_str_out = response.stdout # str:{"jobId":"bd7ac2524f2283e3dd2fa0fa614851ec770a0e47"}
    jobid = re.findall(r"\{\"jobId\"\:\"(.+)\"\}", dict_str_out)[0]
    return jobid

def run_batches(all_uniprotids): # repartition; assigned numbers DON'T map to chunks
    jobid_partition = {}
    partition_number = 1
    for i in range(0, len(all_uniprotids), 50000):
        if i + 50000 < len(all_uniprotids):
            jobid = run_batch(all_uniprotids[i:i+50000]) 
        else:
            jobid = run_batch(all_uniprotids[i:]) 
        jobid_partition[jobid] = partition_number
        partition_number += 1
    print(jobid_partition)
    return jobid_partition

# with open("Data/uniprot_ids.tsv", "r") as f:
#     all_ids_string = f.readlines()[0].strip()
#     uniprot_ids_list = all_ids_string.split("\t")

# job_partition_dict = run_batches(uniprot_ids_list)
# # save to Data/jobid_partition.dict.txt
# with open("Data/jobid_partition.txt", "a") as outfile:
#      print(str(job_partition_dict), file=outfile)

def check_jobstatus(jobid):
    command = "curl -i 'https://rest.uniprot.org/idmapping/status/{}'".format(jobid)
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
        shell = True)
    print(response.stdout)

check_jobstatus("7fb70f75faccae3fdd32f2da6ad13e1076b9973d")
asdf

# retrieve job results
def download_tsvs(jobid, partition):
    url = 'https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{}?compressed=true&fields=accession&format=tsv'.format(jobid) 
    "--output" 
    destination = "Data/Uniprot2NCBI/partition{}.{}.gz".format(str(partition), jobid)
    command = ["curl", url, "--output", destination]
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True)

with open("Data/jobid_partition.dict.txt", "r") as f:
    line = f.readlines()[0].strip()
    job_partition_dict = ast.literal_eval(line)

for jobid, partition in job_partition_dict.items():
    download_tsvs(jobid, partition)