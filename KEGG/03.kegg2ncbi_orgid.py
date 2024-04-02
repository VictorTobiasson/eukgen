import subprocess
import regex as re

#datafile = "Data/taxid_kegg_mappings.br08610.keg"
infile = "Data/kegg_orgs.txt"
with open(infile, "r") as f:
    lines = f.readlines()
    orgs = [line.strip() for line in lines]

baseurl = "https://www.genome.jp/kegg-bin/show_organism?org="
orgs2taxid_dict = {}

mappings = "Data/taxid_mappings.csv"
errorlog = "mappings.log"

for org in orgs:
    lowercase_cds = org.lower()
    url = baseurl + lowercase_cds
    try:
        response = subprocess.run(
                    ["curl", url],
                    text=True,
                    stdout=subprocess.PIPE,
                    check=True)
        out_txt = response.stdout
        tax_id = re.findall(r"https:\/\/www.ncbi.nlm.nih.gov\/Taxonomy\/Browser\/wwwtax.cgi\?mode=Info&id=([0-9]+)\"\>", out_txt)[0]
        orgs2taxid_dict[org] = tax_id
        with open(mappings, "a") as outfile:
            print(",".join([org,tax_id]), file = outfile)
    except:
        with open(errorlog, "a") as logfile:
            print(org, file = logfile)

# output: <tr><td nowrap valign="top"><b>Taxonomy</b></td><td>TAX: 
# <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=36329">36329</a></td></tr>

