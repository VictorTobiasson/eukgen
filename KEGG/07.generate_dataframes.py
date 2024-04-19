import os
import pandas as pd
import numpy as np
import ete3

baseurl_stats = "/data/luojaa/uid_stats/"
baseurl_out = "/data/luojaa/uid_stats_test/"
genes2uid = baseurl_stats + "kegg_genes.mappings.csv" # scraped from KEGG and UIDs mapped using uniprot's idmapping API
uid2kogs = baseurl_stats + "uid2kogs.csv" # scraped from KEGG's internal KOG-level list of UIDs
uid2taxid = baseurl_stats + "uid2taxids.combined.tsv"
uid2descriptor = baseurl_stats + "uid2descriptors.tsv"

genes2uid_df = pd.read_csv(genes2uid).rename(columns={"UNIPROT_ID":"UID"}).set_index("UID")
uid2kogs_df = pd.read_csv(uid2kogs).iloc[:, 1:].rename(columns = {"0":"UID","1":"KOG"}).set_index("UID")
uid2taxid_df = pd.read_csv(uid2taxid, sep = "\t", header = None).rename(columns={0:"UID",1:"TAXID"}).set_index("UID")
uid2desc_df = pd.read_csv(uid2descriptor, sep = "\t", header = None).rename(columns={0:"UID",1:"NAME"}).set_index("UID")

kogs_df = pd.merge(uid2kogs_df, genes2uid_df, on="UID", how = "outer")
kogs_df = pd.merge(kogs_df, uid2taxid_df, on = "UID", how = "outer")
kogs_df = pd.merge(kogs_df, uid2desc_df, on="UID", how = "outer")
kogs_df["KOGID"] = kogs_df["KOG"].fillna(kogs_df["ENTRY"]) # ignore "ENTRY" and "KOG" labels

# get number of missing UIDs per KOG
kogs_df["isnull"] = kogs_df.index.isnull()
null_count = kogs_df.reset_index().loc[:,["isnull", "KOGID"]].groupby("KOGID").sum().astype(int).reset_index().rename(columns={"isnull":"UIDS_MISSING"})
kogs_df = kogs_df.reset_index()
kogs_df = pd.merge(kogs_df, null_count, on = "KOGID", how = "outer")

# label eukaryotic genes
euk_uid2kogs = "/data/luojaa/uid_stats/euk_uid2kogs.csv" # got eukaryotic UIDs from KOG fasta files (forgot to record when fastas were create)
euk_uid2kogs_df = pd.read_csv(euk_uid2kogs, header = None).rename(columns = {0:"UID", 1:"KOGID"}).set_index("UID")
euk_uid2kogs_df["ISEUK"] = [True] * len(euk_uid2kogs_df)
euk_uid2kogs_df = euk_uid2kogs_df.loc[:,"ISEUK"]
kogs_df_tmp = kogs_df.set_index("UID")
kogs_df_iseuk = pd.merge(kogs_df_tmp, euk_uid2kogs_df, on = "UID", how="outer")

kogs_df_iseuk.fillna({"ISEUK":False}, inplace = True)
euk_count = kogs_df_iseuk.reset_index().loc[:,["ISEUK", "KOGID"]].groupby("KOGID").sum().astype(int).reset_index().rename(columns={"ISEUK":"EUKCOUNT"}) # count the eukaryotic genes per KOG
kogs_df_iseuk = kogs_df_iseuk.reset_index()
kogs_df_summary = pd.merge(kogs_df_iseuk, euk_count, on = "KOGID", how = "outer") # intermediate file

# BRANCH 1: map taxids to taxonomic names
ncbi = ete3.NCBITaxa()
taxids = list(kogs_df_summary[~kogs_df_summary["TAXID"].isnull()]["TAXID"])
taxid2name = ncbi.get_taxid_translator(taxids)
def translate_taxid(name, dic):
    try:
        return dic[name]
    except:
        return np.nan
kogs_df_taxa = kogs_df_summary["TAXID"].apply(lambda x: translate_taxid(x, taxid2name))

kogs_df_crucial = kogs_df_summary.loc[:, ["UID", "KOGID", "NAME", "ALIAS", "KEGG_ORG", "ISEUK"]]
kogs_df_crucial["TAXID"] = kogs_df_summary["TAXID"].fillna(0).astype(int)
kogs_df_crucial["SPECIES"] = kogs_df_taxa
kogs_df_searchable = kogs_df_crucial.fillna("none")

# BRANCH 2: get kog (or fasta) level stats
fasta_stats_df = kogs_df_summary.loc[:,["KOGID", "UID", "UIDS_MISSING", "EUKCOUNT"]].groupby(["KOGID","UIDS_MISSING", "EUKCOUNT"]).count().reset_index().rename(columns={"UID":"UID_COUNT"}).set_index("KOGID")
cluster_size = pd.DataFrame(kogs_df_summary.loc[:,["KOGID"]].groupby("KOGID").size()).rename(columns={0:"KOGSIZE"})
fasta_stats_out = pd.merge(cluster_size, fasta_stats_df, on = "KOGID", how = "outer").reset_index()
fasta_stats_out["EUK_FRACTION"] = 100 * fasta_stats_out["EUKCOUNT"]/fasta_stats_out["UID_COUNT"]

# format "CATEGORIES" into a stack
pathways_path = "/data/luojaa/kegg/kegg_pathways.58345.tsv" # from KEGG scraping
reactions_path = "/data/luojaa/kegg/kegg_reactions.tsv"
modules_path = "/data/luojaa/kegg/kegg_modules.tsv"
kegg_pathways = pd.read_csv(pathways_path, sep = "\t").loc[:,["ENTRY", "PATHWAY_ID", "PATHWAY_NAME"]].rename(
    columns = {"PATHWAY_ID":"CATEGORY_ID", "PATHWAY_NAME":"CATEGORY_NAME", "ENTRY":"KOGID"})
kegg_rxns = pd.read_csv(reactions_path, sep = "\t").rename(columns = {"REACTION_ID":"CATEGORY_ID", "REACTION_NAME":"CATEGORY_NAME", "ENTRY":"KOGID"})
kegg_modules = pd.read_csv(modules_path, sep = "\t").rename(columns = {"MODULE_ID":"CATEGORY_ID", "MODULE_NAME":"CATEGORY_NAME", "ENTRY":"KOGID"})
kegg_categories_df = pd.concat([kegg_pathways, kegg_rxns, kegg_modules], axis = 0).dropna().reset_index().iloc[:, 1:] # drop misread index column

# elongate categories stack for searchability
kegg_categories = kegg_categories_df.copy(deep=True)
kegg_categories["KOGIDS"] = kegg_categories.groupby(["CATEGORY_NAME", "CATEGORY_ID"])["KOGID"].transform(lambda x: ",".join(x))
kegg_categories["NUM_KOGS"] = kegg_categories.groupby(["CATEGORY_ID"]).transform(lambda x: len(x))["KOGID"]
kegg_categories_groups = kegg_categories.loc[:,["CATEGORY_NAME", "CATEGORY_ID", "KOGIDS", "NUM_KOGS"]].drop_duplicates().dropna()

# get kogs with proteins in both euk and prok; NOTE: some of these proteins don't map to our sequence db (yet)
uid2taxid_iseuk = base_stats + "uid2taxids2iseuk.combined.tsv"
uid2taxid_iseuk_df = pd.read_csv(uid2taxid_iseuk, sep = "\t", header = None).rename(columns={0:"UID",1:"TAXID",2:"iseuk"}).set_index("UID")
kogs_domains = pd.merge(uid2taxid_iseuk_df, uid2kogs_df, on = "UID", how = "left")
kogs_domains["isprok"] = ~kogs_domains["iseuk"] 
has_prok = kogs_domains.reset_index().loc[:,["KOG","isprok"]].groupby("KOG").apply(lambda x: x.eq(True).any()).loc[:,["isprok"]].rename(columns={"isprok":"hasprok"})
has_euk = kogs_domains.reset_index().loc[:,["KOG","iseuk"]].groupby("KOG").apply(lambda x: x.eq(True).any()).loc[:,["iseuk"]].rename(columns={"iseuk":"haseuk"})
kogs_domains = kogs_domains.reset_index().set_index("KOG")
kogs_hasdomains = pd.merge(kogs_domains, has_prok, on="KOG", how = "left")
kogs_hasdomains = pd.merge(kogs_hasdomains, has_euk, on="KOG", how = "left")
kogs_hasboth = kogs_hasdomains.loc[:,["hasprok", "haseuk"]].reset_index().drop_duplicates()
kogs_hasboth["hasboth"] = kogs_hasboth["hasprok"] & kogs_hasboth["haseuk"]
kogs_hasboth["haseither"] = kogs_hasboth["hasprok"] | kogs_hasboth["haseuk"]
KOGS_w_both = kogs_hasboth[kogs_hasboth["hasboth"] == True]["KOG"]

# OUTPUTS
kogs_df.to_csv(baseurl_out + "kogs_df.tsv", sep = "\t", index=None) # intermediate file. no euk info. 
kogs_df_crucial.to_csv(baseurl_out + "kogs_df_crucial.tsv", sep = "\t", index=None) # includes iseuk and species, excluedes eukcount and nullcount
kogs_df_searchable.to_csv(baseurl_out + "kogs_df_searchable.tsv", sep = "\t", index=None) # nans filled in
fasta_stats_out.to_csv(baseurl_out + "cluster_stats.csv", sep = ",", index = None) 
kegg_categories_df.to_csv(baseurl_out + "kegg_categories.tsv", sep = "\t", index = None) 
kegg_categories_groups.to_csv(baseurl_out + "kegg_categories_searchable.tsv", sep = "\t", index = None)
KOGS_w_both.to_csv(baseurl_out + "KOGs_w_both.csv")