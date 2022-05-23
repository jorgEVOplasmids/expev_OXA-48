#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 12:59:02 2022

@author: jorge

This script takes the output of mut_screening.py and renames genes annotated
with PGAP ids (pgaptmp_number) changing it by the breseq output description.

"""

import pandas as pd

strain = "K209"

all_muts = pd.read_excel("/home/jorge/Documents/TFM/clinical_samples/R_analysis/all_mutations_heatmaps/"+strain+"_all_muts.xlsx")
breseq_output_muts = pd.read_excel("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/summary_table_"+strain+"_p.xlsx", sheet_name = "mutations")
breseq_output_NJ = pd.read_excel("/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs/summary_table_"+strain+"_p.xlsx", sheet_name = "JC full")

all_muts.index = all_muts["Unnamed: 0"]
all_muts = all_muts.drop("Unnamed: 0", axis = 1)

rownames = list(all_muts.index.values)

id_genename_dic = {}

for index, row in breseq_output_muts.iterrows():
    if "pgaptmp" in row[5] and row [5] not in id_genename_dic.keys():
        id_genename_dic[row[5]] = row[6] + "*" # the * marks the gene as mutated by a SNP or Indel
        
for index, row in breseq_output_NJ.iterrows():
    if "pgaptmp" in row[10] and row [10] not in id_genename_dic.keys():
        id_genename_dic[row[10]] = row[11]
    if "pgaptmp" in row[13] and row [13] not in id_genename_dic.keys():
        id_genename_dic[row[13]] = row[14]

renamed_all_muts = all_muts

for i, row in renamed_all_muts.iterrows():
    if "pgaptmp" in i:
        renamed_all_muts = renamed_all_muts.rename(index = {i: id_genename_dic[i]})
        
# drop hypothetical proteins to avoid problems later with duplicated rownames
"""
for i, row in renamed_all_muts.iterrows():
    if "hypothetical protein" == i or "hypothetical protein/hypothetical protein" == i:
        print(i)
        renamed_all_muts = renamed_all_muts.drop(labels = i)
    else:
        renamed_all_muts = renamed_all_muts
"""
renamed_all_muts.to_excel("/home/jorge/Documents/TFM/clinical_samples/R_analysis/all_mutations_heatmaps_renamed/"+strain+"_all_muts_renamed.xlsx")
