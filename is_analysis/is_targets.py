#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:48:16 2022

@author: jorge

This script takes a breseq excel output parsed with the breseq_index_parser.py script
and gets genes targeted by NJ evidences due to pOXA-48 IS movements for further analysis.

"""

import os
import pandas as pd

DirName = "/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs"
files = os.listdir(DirName)
strain = "K209"
list_count_genes = []
oxa_is = ["IS1 family transposase/carbapenem‑hydrolyzing class D beta‑lactamase OXA‑48",
          "IS1 family transposase/DGQHR domain‑containing protein",
          "DGQHR domain‑containing protein",
          "IS4‑like element IS10A family transposase",
          "DGQHR domain‑containing protein/IS1 family transposase"]

for file in files:
    file_whole = DirName + "/" + file
    if strain in file: # for test
        xlsx = pd.ExcelFile(file_whole)
        muts = pd.read_excel(xlsx, "mutations")
        NJ = pd.read_excel(xlsx, "JC full")
        # For the NJ sheet we have to do some extra tasks
        # Define the necessary variables
        pos_to_drop = []
        NJ_contig = NJ["Seq ID 1"]
        NJ_position = NJ["Position 1"]
        NJ_gene1 = NJ["Gene 1"]
        NJ_gene2 = NJ["Gene 2"]
        NJ_prod1 = NJ["Product 1"]
        NJ_prod2 = NJ["Product 2"]
        NJ_pOXA_Anc = NJ.iloc[:, 15]
        NJ_pOXA_1 = NJ.iloc[:, 16]
        NJ_pOXA_2 = NJ.iloc[:, 17]
        NJ_pOXA_3 = NJ.iloc[:, 18]
        NJ_nopOXA_Anc = NJ.iloc[:, 19]
        NJ_nopOXA_1 = NJ.iloc[:, 20]
        NJ_nopOXA_2 = NJ.iloc[:, 21]
        NJ_nopOXA_3 = NJ.iloc[:, 22]
        NJ_pOXA_Ab_1 = NJ.iloc[:, 23]
        NJ_pOXA_Ab_2 = NJ.iloc[:, 24]
        NJ_pOXA_Ab_3 = NJ.iloc[:, 25]
        genes_pOXA_1 = {}
        genes_pOXA_2 = {}
        genes_pOXA_3 = {}
        genes_pOXA_Ab_1 = {}
        genes_pOXA_Ab_2 = {}
        genes_pOXA_Ab_3 = {}
        # Firstly, get a list of the targeted positions in the chromosome in the ancestors
        for i in range(len(NJ)):
            if NJ_pOXA_Anc[i] != 0 or NJ_nopOXA_Anc[i] != 0:
                pos_to_drop.append(NJ_position[i])
        # Now, get the genes targeted by the new junctions of pOXA-48 IS
        # If the gene is targeted by both ends of the NJ it is only counted once
        for i in range(len(NJ)):
            if NJ_pOXA_1[i] != 0 and NJ_pOXA_1[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if str(NJ_prod2[i]) in oxa_is:
                        if NJ_gene1[i] not in genes_pOXA_1:
                            genes_pOXA_1[NJ_gene1[i]] = NJ_pOXA_1[i]
                        else:
                            genes_pOXA_1[NJ_gene1[i]] += NJ_pOXA_1[i]
                    if  str(NJ_prod1[i]) in oxa_is:
                        if NJ_gene2[i] not in genes_pOXA_1:
                            genes_pOXA_1[NJ_gene2[i]] = NJ_pOXA_1[i]
                        else:
                            genes_pOXA_1[NJ_gene2[i]] += NJ_pOXA_1[i]
            if NJ_pOXA_2[i] != 0 and NJ_pOXA_2[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if str(NJ_prod2[i]) in oxa_is:
                        if NJ_gene1[i] not in genes_pOXA_2:
                            genes_pOXA_2[NJ_gene1[i]] = NJ_pOXA_2[i]
                        else:
                            genes_pOXA_2[NJ_gene1[i]] += NJ_pOXA_2[i]
                    if  str(NJ_prod1[i]) in oxa_is:
                        if NJ_gene2[i] not in genes_pOXA_2:
                            genes_pOXA_2[NJ_gene2[i]] = NJ_pOXA_2[i]
                        else:
                            genes_pOXA_2[NJ_gene2[i]] += NJ_pOXA_2[i]
            if NJ_pOXA_3[i] != 0 and NJ_pOXA_3[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if str(NJ_prod2[i]) in oxa_is:
                        if NJ_gene1[i] not in genes_pOXA_3:
                            genes_pOXA_3[NJ_gene1[i]] = NJ_pOXA_3[i]
                        else:
                            genes_pOXA_3[NJ_gene1[i]] += NJ_pOXA_3[i]
                    if  str(NJ_prod1[i]) in oxa_is:
                        if NJ_gene2[i] not in genes_pOXA_3:
                            genes_pOXA_3[NJ_gene2[i]] = NJ_pOXA_3[i]
                        else:
                            genes_pOXA_3[NJ_gene2[i]] += NJ_pOXA_3[i]
            if NJ_pOXA_Ab_1[i] != 0 and NJ_pOXA_Ab_1[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if str(NJ_prod2[i]) in oxa_is:
                        if NJ_gene1[i] not in genes_pOXA_Ab_1:
                            genes_pOXA_Ab_1[NJ_gene1[i]] = NJ_pOXA_Ab_1[i]
                        else:
                            genes_pOXA_Ab_1[NJ_gene1[i]] += NJ_pOXA_Ab_1[i]
                    if  str(NJ_prod1[i]) in oxa_is:
                        if NJ_gene2[i] not in genes_pOXA_Ab_1:
                            genes_pOXA_Ab_1[NJ_gene2[i]] = NJ_pOXA_Ab_1[i]
                        else:
                            genes_pOXA_Ab_1[NJ_gene2[i]] += NJ_pOXA_Ab_1[i]
            if NJ_pOXA_Ab_2[i] != 0 and NJ_pOXA_Ab_2[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if str(NJ_prod2[i]) in oxa_is:
                        if NJ_gene1[i] not in genes_pOXA_Ab_2:
                            genes_pOXA_Ab_2[NJ_gene1[i]] = NJ_pOXA_Ab_2[i]
                        else:
                            genes_pOXA_Ab_2[NJ_gene1[i]] += NJ_pOXA_Ab_2[i]
                    if  str(NJ_prod1[i]) in oxa_is:
                        if NJ_gene2[i] not in genes_pOXA_Ab_2:
                            genes_pOXA_Ab_2[NJ_gene2[i]] = NJ_pOXA_Ab_2[i]
                        else:
                            genes_pOXA_Ab_2[NJ_gene2[i]] += NJ_pOXA_Ab_2[i]
            if NJ_pOXA_Ab_3[i] != 0 and NJ_pOXA_Ab_3[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if str(NJ_prod2[i]) in oxa_is:
                        if NJ_gene1[i] not in genes_pOXA_Ab_3:
                            genes_pOXA_Ab_3[NJ_gene1[i]] = NJ_pOXA_Ab_3[i]
                        else:
                            genes_pOXA_Ab_3[NJ_gene1[i]] += NJ_pOXA_Ab_3[i]
                    if  str(NJ_prod1[i]) in oxa_is:
                        if NJ_gene2[i] not in genes_pOXA_Ab_3:
                            genes_pOXA_Ab_3[NJ_gene2[i]] = NJ_pOXA_Ab_3[i]
                        else:
                            genes_pOXA_Ab_3[NJ_gene2[i]] += NJ_pOXA_Ab_3[i]
        # And finally include 
        list_count_genes.append(genes_pOXA_1)
        list_count_genes.append(genes_pOXA_2)
        list_count_genes.append(genes_pOXA_3)
        list_count_genes.append(genes_pOXA_Ab_1)
        list_count_genes.append(genes_pOXA_Ab_2)
        list_count_genes.append(genes_pOXA_Ab_3)

muts = pd.DataFrame(list_count_genes)
muts = muts.T
muts = muts.fillna(0)
muts = muts.rename(columns = ({0:"pOXA-48 1", 1:"pOXA-48 2", 2: "pOXA-48 3",
                               3: "pOXA-48+Ab 1", 4: "pOXA-48+Ab 2", 5: "pOXA-48+Ab 3"}))

muts.to_excel("/home/jorge/Documents/TFM/clinical_samples/R_analysis/oxa_is_targets_xlsx/"+str(strain)+"_oxa_is_targets_xlsx.xlsx")
