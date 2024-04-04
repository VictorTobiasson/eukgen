import subprocess
import regex as re
import sys
from collections import defaultdict

#knumbers_unfiltered = ["K" + str((5-len(str(i)))*"0") + str(i) for i in range(1000)]

def concat_clusters(group_no, knumbers):
    knumbers_range = knumbers[(group_no - 1) * 100: group_no * 100]
    paths = ["data/kog_uid/{}.tsv".format(knumber) for knumber in knumbers_range]
    outfilename = "data/uid_groups/group{}.tsv".format(group_no)
    with open(outfilename, "a") as outfile:
        for path in paths:
            try:
                with open(path, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        uid = line.strip().split("\t")[1]
                        print(uid, file=outfile)
            except:
                continue
def assign_batches(group_no):
    batch_ids = defaultdict(list)
    batch_id, ind = 1, 1
    infilename = "data/uid_groups/group{}.tsv".format(group_no)
    with open(infilename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if ind % 10 == 0:
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
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
        shell=True)
    print(command)
    dict_str_out = response.stdout # str:{"jobId":"62910de2cb548388eac8e873394b006e9c2c16de"}
    print(dict_str_out)
    jobid = re.findall(r"\{\"jobId\"\:\"(.+)\"\}", dict_str_out)[0]
    return jobid

def run_batches(groups, section):
    for group_no in groups:
        concat_clusters(group_no, knumbers_unfiltered)
        batches = assign_batches(group_no)
        with open("data/batch_jobids/section{}.tsv".format(section), "a") as outfile:
            for batch in batches.keys():
                uids = batches[batch]
                print(len(uids))
                jobid = run_batch(uids)
                batch_id = "{}.{}".format(group_no, batch)
                print("\t".join([batch_id, jobid]), file=outfile)
run_batches([1, 2], "testsection")



