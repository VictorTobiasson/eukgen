import subprocess
import regex as re
import sys
from collections import defaultdict

knumbers_unfiltered = ["K" + str((5-len(str(i)))*"0") + str(i) for i in range(27110)]

start_group, end_group, section = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

def concat_clusters(group_no, knumbers):
    knumbers_range = knumbers[(group_no - 1) * 100: group_no * 100]
    paths = ["/data/luojaa/kegg/kog_uid/{}.tsv".format(knumber) for knumber in knumbers_range]
    outfilename = "/data/luojaa/kegg/uid_groups/group{}.tsv".format(group_no)
    logfilename = "/data/luojaa/kegg/uidsdisordered.log"
    with open(outfilename, "a") as outfile, open(logfilename, "a") as logfile:
        for path in paths:
            try:
                with open(path, "r") as f:
                    uidprev = "00000"
                    lines = f.readlines()
                    for line in lines:
                        uid = line.strip().split("\t")[1]
                        if uid < uidprev:
                            print(group_no, file=logfile)
                        print(uid, file=outfile)
                        uid = uidprev
            except:
                continue
def assign_batches(group_no):
    batch_ids = defaultdict(list)
    batch_id, ind = 1, 1
    infilename = "/data/luojaa/kegg/uid_groups/group{}.tsv".format(group_no)
    with open(infilename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if ind % 10000 == 0:
                batch_id += 1
            batch_ids[batch_id] += [line.strip()]
            ind += 1
    return batch_ids
def run_batch(ids):
    requesturl = "https://rest.uniprot.org/idmapping/run"
    from_db = "\'from=\"UniProtKB_AC-ID\"\'"
    to_db = "\'to=\"UniRef100\"\'"
    idstr = "\'ids=\"" + ",".join(ids) + "\"\'"
    command = "curl --request POST {} --form {} --form {} --form {}".format(requesturl, from_db, to_db, idstr)
    response = subprocess.run(command, text=True, stdout=subprocess.PIPE, check=True, shell=True)
    print(command)
    dict_str_out = response.stdout # str:{"jobId":"62910de2cb548388eac8e873394b006e9c2c16de"}
    print(dict_str_out)
    jobid = re.findall(r"\{\"jobId\"\:\"(.+)\"\}", dict_str_out)[0]
    return jobid

def run_batches(groups, section):
    for group_no in groups:
        concat_clusters(group_no, knumbers_unfiltered)
        batches = assign_batches(group_no)
        with open("/data/luojaa/kegg/batch_jobids/section{}.tsv".format(section), "a") as outfile:
            for batch in batches.keys():
                uids = batches[batch]
                print(len(uids))
                jobid = run_batch(uids)
                batch_id = "{}.{}".format(group_no, batch)
                print("\t".join([batch_id, jobid]), file=outfile)
groups = [i for i in range(start_group, end_group)]
run_batches(groups, section)



