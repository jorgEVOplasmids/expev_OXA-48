#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:23:01 2022

@author: jorgesd

### Script to count the IS movements in each experimental group. It summarizes,
from the breseq xlsx parsed output, the number of IS rearrangements produced 
with/without pOXA, and with/without antibiotic ###

"""

import os
import pandas as pd


DirName = "/home/jorgesd/Documents/TFM/R_analysis/xlsx_breseqs"

files = os.listdir(DirName)

summary_list = []

for file in files:
    strain = file.split("_")[2]
    file_whole = DirName + "/" + file
    xlsx = pd.ExcelFile(file_whole)
    NJ = pd.read_excel(xlsx, "JC full")
    pos_to_drop = []
    strain_stats = []
    NJ_contig1 = NJ["Seq ID 1"]
    NJ_contig2 = NJ["Seq ID 2"]
    NJ_position = NJ["Position 1"]
    NJ_annot1 = NJ["Product 1"]
    NJ_annot2 = NJ["Product 2"]
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
    pOXA_IS = 0 # will contain the pOXA without Ab rearrangements
    nopOXA_IS = 0 # will contain the no pOXA rearrangements
    pOXA_Ab_IS = 0 # will contain the pOXA+Ab rearrangements
    OXA_Ab_rear = 0
    OXA_rear = 0
    # Get a list of the positions targeted by the ancestors samples to remove afterwards
    for i in range(len(NJ)):
        if NJ_contig1[i] == 1 or NJ_contig1[i] == "contig_1":
            if NJ_pOXA_Anc[i] != 0 or NJ_nopOXA_Anc[i] != 0 or NJ_nopOXA_1[i] != 0 or NJ_nopOXA_2[i] != 0 or NJ_nopOXA_3[i] != 0:
                pos_to_drop.append(NJ_position[i])
    # Now, count the IS movements by condition
    for i in range(len(NJ)):
        if NJ_contig1[i] == 1 or NJ_contig1[i] == "contig_1":
            if NJ_pOXA_1[i] != 0 or NJ_pOXA_2[i] != 0 or NJ_pOXA_3[i] != 0 and NJ_position[i] not in pos_to_drop: #for samples with OXA and no Ab
                if "IS" in NJ_annot1[i] or "IS" in NJ_annot2[i] or "DGQHR" in NJ_annot1[i] or "DGQHR" in NJ_annot2[i]:
                    pOXA_IS += 1
                if strain == "CF12" or strain == "CF13":
                    if NJ_contig2[i] == 2:
                        OXA_rear += 1
                if strain == "H53" or strain == "K091" or strain == "K147" or strain == "K163":
                    if NJ_contig2[i] == 3:
                        OXA_rear += 1
                if strain == "K153":
                    if NJ_contig2[i] == "contig_2":
                        OXA_rear += 1
                if strain == "K209":
                    if NJ_contig2[i] == 4:
                        OXA_rear += 1
            if NJ_nopOXA_1[i] != 0 or NJ_nopOXA_2[i] != 0 or NJ_nopOXA_3[i] != 0 and NJ_position[i] not in pos_to_drop: #for samples without OXA
                if "IS" in NJ_annot1[i] or "IS" in NJ_annot2[i] or "DGQHR" in NJ_annot1[i] or "DGQHR" in NJ_annot2[i]:
                    nopOXA_IS += 1
            if NJ_pOXA_Ab_1[i] != 0 or NJ_pOXA_Ab_2[i] != 0 or NJ_pOXA_Ab_3[i] != 0 and NJ_position[i] not in pos_to_drop: #for samples without OXA
                if "IS" in NJ_annot1[i] or "IS" in NJ_annot2[i] or "DGQHR" in NJ_annot1[i] or "DGQHR" in NJ_annot2[i]:
                    pOXA_Ab_IS += 1
                if strain == "CF12" or strain == "CF13":
                    if NJ_contig2[i] == 2:
                        OXA_Ab_rear += 1
                if strain == "H53" or strain == "K091" or strain == "K147" or strain == "K163":
                    if NJ_contig2[i] == 3:
                        OXA_Ab_rear += 1
                if strain == "K153":
                    if NJ_contig2[i] == "contig_2":
                        OXA_Ab_rear += 1
                if strain == "K209":
                    if NJ_contig2[i] == 4:
                        OXA_Ab_rear += 1
    # print(strain, nopOXA_IS, pOXA_IS, pOXA_Ab_IS, OXA_rear, OXA_Ab_rear)
    strain_stats.append(strain)
    strain_stats.append(nopOXA_IS)
    strain_stats.append(pOXA_IS)
    strain_stats.append(pOXA_Ab_IS)
    strain_stats.append(pOXA_IS+pOXA_Ab_IS)
    strain_stats.append(OXA_rear)
    strain_stats.append(OXA_Ab_rear)
    summary_list.append(strain_stats)
    
# print(summary_list)
final_df = pd.DataFrame(summary_list)
"""column_names = ["Strain", "IS rearrangements (no pOXA-48)", "IS rearrangements (pOXA-48 no Ab)",
                "IS rearrangements (pOXA-48+Ab)", "IS rearrangements (total pOXA-48)",
                "IS rearrangements due to pOXA-48", "IS rearrangements due to pOXA-48 in presence of Ab"]"""
# final_df = final_df.reindex(columns = column_names)
final_df = final_df.reindex([1,0,4,5,6,2,7,3])
final_df = final_df.T
# print(final_df)

final_df.to_excel("IS_rearrangements.xlsx")