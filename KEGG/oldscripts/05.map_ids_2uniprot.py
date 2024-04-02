import subprocess
import regex as re
import ast

# map chunks to kegg ids
def read_ids():
    chunk_keggids = {}
    ## get ids from chunks
    chunk_directory = "Data/KEGG_CDS_chunks/"
    chunk_paths = [chunk_directory + "chunk{}.txt".format(str(i)) for i in range(72, 272)] #did 71 chunks manually
    for chunk in chunk_paths:
        chunk_number = int(re.findall(r"chunk([0-9]+).txt", chunk)[0])
        keggids = []
        with open(chunk, "r") as f:
            lines = f.readlines()
            for line in lines: # chunks should have one line
                some_ids = line.strip().split(" ")
                keggids += some_ids
        chunk_keggids[chunk_number] = keggids
    return chunk_keggids

#chunk_keggids = read_ids()


## run batch requests on ids
#ids = ["hsa:55615","hsa:79899"]
def run_batch(ids):
    requesturl = "https://rest.uniprot.org/idmapping/run"
    from_db = "\'from=\"KEGG\"\'"
    to_db = "\'to=\"UniProtKB\"\'"
    idstr = "\'ids=\"" + ",".join(ids) + "\"\'"
    command = "curl --form {} --form {} --form {} {}".format(from_db, to_db, idstr, requesturl)
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
        shell = True)
    dict_str_out = response.stdout # str:{"jobId":"62910de2cb548388eac8e873394b006e9c2c16de"}
    jobid = re.findall(r"\{\"jobId\"\:\"(.+)\"\}", dict_str_out)[0]
    return jobid

def run_batches(chunk_keggids):
    jobid_chunk = {}
    for chunk_number in list(chunk_keggids.keys()):
        keggids = chunk_keggids[chunk_number]
        jobid_1 = run_batch(keggids[:50000]) #split batch in half to overcome arg_max bash limits
        jobid_2 = run_batch(keggids[50000:])
        jobid_chunk[jobid_1] = chunk_number
        
        jobid_chunk[jobid_2] = chunk_number
    return jobid_chunk # saved output to Data/jobid_chunk.dict.txt

def check_jobstatus(jobid):
    command = "curl -i 'https://rest.uniprot.org/idmapping/status/{}'".format(jobid)
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
        shell = True)
    print(response.stdout)

#check_jobstatus("387ca94b78a9f27d00420620eab1d6d4d1b83a2b")

# retrieve job results
def download_tsvs(jobid, chunk):
    url = 'https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{}?compressed=true&fields=accession&format=tsv'.format(jobid) 
    "--output" 
    destination = "Data/Kegg2Uniprot/chunk{}.{}.gz".format(str(chunk), jobid)
    command = ["curl", url, "--output", destination]
    response = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        check=True)

with open("Data/jobid_chunk.dict.txt", "r") as f:
    lines = f.readlines()
    jobid_dict = ast.literal_eval(lines[0].strip())
    print(jobid_dict)

for jobid, chunk in jobid_dict.items():
    download_tsvs(jobid, chunk)


# get output from job
# job_ids = []
# job_ids = ["6baeb09ef587b63780e4abe18e7be01a5608dbef", 
#            "4efd8e3e581ed6a31597ae2c54be04e97ad2e971", 
#            "9841a8289177920d9c4ca0eb04c018a0a59bcc9d", 
#            "2412d8fae754eb7fe64b9674d56d13ad72ab8189",
#            "c9afce5181131c49674a6e7c36e27e120f3c5ddf",
#            "8f57f0cb597ead8f898ff16eb5f01c5f48426c73",
#            "797a99f3f3212b601564744c38ae4da4321a693a",
#            "18836929be4b22a17668704a6486f7ffe0000a01",
#            "5cb8b8f48fab688bbe790ec6be9f55e5ac68b112",
#            "3b1739cafb68193978342fa7b38c2644ccdcc4d4",
#            "5b13bd783c7edf62cadba1b7d8fa2fa5eefcee26",
#            "83ad1bccdd3244a6243044a557e95f4ef1c48119",
#            "115a8abd77e410acb1baf383cd06efdbbceed79e",
#            "5d8130ed701ac2e18fce19da6b82b83637c540b4",
#            "60d51a7608784444f2b6d16ad4c6b1d873df2d59",
#            "036be1a72a31b9b291c94ad3f378eeec9affc4fd",
#            "716f6fd904577f03a16bc44308217dba1b24a6e1",
#            "aeb9a903348ad044c73d0907b32b3379dec86696",
#            "0cb706db96d291e5e31fc7398ded95df3d8002aa",
#            "7892c33360ec792b1da8ea46f9599c767e2e88ee",
#            "1bc65864bcdf703e6754d29cf73feb927ac0e45a",
#            "2c19bc2793a96feab43701b7ed7e873e0e4e2a0f",
#            "7819f492dbf97e43609fee7d6c3442dd78fdfd1a",
#            "5829ffd7b42b76329800c91c1c76564b5cfa5ba0",
#            "729ecf645555c85c2798f315a1f9bc960c656ecc",
#            "a4f51e9a455e315b43170794cd3bf1ed62a1f02e",
#            "808433bcd52cab5d6ed4a4b2efc31bed6c977096",
#            "f3c02b13f5e93b960e8ef349d6806e03780a03fb",
#            "5d81e2b97dac26f4eebc3fe454cedc423facf6f7",
#            "1d5d261803f2e0de35a26d3d63068d0119b128bd",
#            "ba442195748bf1d08caf74e61688a120d725e658",
#            "59a08f20624ed9e5084d7f7e9bccea1802933a4f",
#            "0f32f54562e1956c851da60d731866ea86e0dbdd",
#            "1bebe971fb2134e779427fed7275368ae3b1b2da",
#            "45456881a595504f8633ffd73bb541b009356816",
#            "bd0c7ec48aec2c4aa36019e852f0d52a7390af91",
#            "c6243b5acb077dbb032b8e956ed08b0fb6620bfb",
#            "30408731142f6fb1b002c6152f92da727999d05e",
#            "fa97d2d1d9e2b449616a7eca88e7524e82c7a1cd",
#            "b1f15a077c39258155f2d96861cf4230dea78562",
#            "4ed60ec640e1dd0805c4d5b6268b367cbf25e473",
#            "dadf8b99424b84b322cc5428329edc167749f4d2",
#            "fed5cbfd1aac24de1ee3673ef88810c3ce7b00d5",
#            "f2e6ddf3c96b1ba5950e74a2d58186d1d611b0ae",
#            "bbb0c72ed6b02883e2da4a4d24d6d9a621d0d12e",
#            "32ca98ccb3cac737ec7ce14395ce927b124c38a7",
#            "ae8ac4add788734c5e381a13b4df032f3c337219",
#            "aac55f580798562bf3870b2845e0afe721d26c03",
#            "f8bdc5f3c7ee8b714138d6c7cfddebcaf2d2bc3a"]


# chunk = 23
# for jobid in job_ids:
#     url = 'https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{}?compressed=true&fields=accession&format=tsv'.format(jobid) 
#     destination = "Data/Kegg2Uniprot/chunk{}.{}.gz".format(str(chunk), jobid)
#     command = ["curl", url, "--output", destination]
#     response = subprocess.run(
#         command,
#         text=True,
#         stdout=subprocess.PIPE,
#         check=True)
#     chunk += 1
