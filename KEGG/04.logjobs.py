import regex as re
import sys
import subprocess

def check_jobstatus(jobid):
    command = f"curl -i 'https://rest.uniprot.org/idmapping/status/{jobid}'"
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
        shell = True)
    return response.stdout

section=sys.argv[1]
jobid_file = f"/data/luojaa/kegg/batch_jobids/section{section}.tsv"
job_log = f"/data/luojaa/kegg/section{section}.log"
fcount, ncount, rcount = 0,0,0
with open(jobid_file, "r") as f, open(job_log, "a") as logfile:
    lines = f.readlines()
    for line in lines:
        batch, jobid = line.strip().split("\t")
        status_out = check_jobstatus(jobid)
        status = re.findall(r"jobStatus\":\"(.*)\"}", status_out)[0]
        print(status)
        if status =="FINISHED":
            fcount += 1
        elif status == "NEW":
            ncount +=1
        elif status == "RUNNING":
            rcount += 1
    print(",".join([count, ncount, rcount]), file=logfile)