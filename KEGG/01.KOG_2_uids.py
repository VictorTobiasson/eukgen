import requests
from bs4 import BeautifulSoup
from pprint import pprint
import regex as re
from multiprocessing import Pool
import sys

base_url = "https://www.kegg.jp/kegg-bin/uniprot_list?ko=" # append KO group to end of url
def transform_int2knumber(knint):
    knstr = str(knint)
    while len(knstr) < 5:
        knstr = "0" + knstr
    return "K" + knstr
not_kog_ids = []
with open("not_knumber.log", "r") as f:
    lines = f.readlines()
    for line in lines:
        knumber = line.strip()
        not_kog_ids.append(knumber)
kog_ids_raw = [transform_int2knumber(k) for k in range(1, 27110)]
kog_ids = [k for k in kog_ids_raw if k not in not_kog_ids]

def extract_uid_htmlrow(rowtext):
    uid = re.findall(r'">(.*)<\/a', rowtext)[0]
    return uid
def extract_descriptor_htmlrow(rowtext, altformat = False):
    if altformat:
        descriptor = re.findall(r'[^\n].*$', descriptor_row)[0].split("\xa0")[-1]
        return descriptor
    descriptor = re.findall(r';\s(.*)', rowtext)[0] 
    return descriptor

def get_uids(kogid): # also gets descriptor for KOG
    uids = []
    url = base_url + kogid
    data = requests.get(url)
    print("completed request for {}".format(kogid))
    html = BeautifulSoup(data.text, 'html.parser')
    try: 
        descriptor_row = html.div.text
    except:
        with open("not_knumber_anymore.log", "a") as outfile:
            print(kogid, file=outfile)
            return
    try:
        descriptor = extract_descriptor_htmlrow(descriptor_row)
    except: # row is formated without ";" separating the gene alias and descriptor
        descriptor = extract_descriptor_htmlrow(descriptor_row, True)
        
    for i in range(2, len(html.find_all("td")), 4):
        rowtext = str(html.find_all("td")[i])
        uid = extract_uid_htmlrow(rowtext)
        uids.append(uid)
    output = "data/kog_uid/{}.tsv".format(kogid)
    with open(output, "a") as outfile:
        for uid in uids:
            print("\t".join([kogid, uid, descriptor]), file=outfile)

section, num_sections = int(sys.argv[1]), int(sys.argv[2])
#print(section)
section_indices = [i for i in range(0, len(kog_ids), len(kog_ids) // num_sections)]
start_range = section_indices[section-1]
if section == num_sections:
    end_range = len(kog_ids)
else:
    end_range = section_indices[section]

start_knumber = int(sys.argv[3])
    
for kog_id in kog_ids[start_range:end_range]:
   # print(kog_id)
    if int(kog_id[1:]) >= start_knumber:
        get_uids(kog_id)

