#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 11:48:28 2022

@author: jorgesd

@title: Mutation Spectrum Script v.4

Final version of the mutation spectrum script. Extracts info for each
sample treated individually and filters out the mutations present in samples
with pOXA-48 and in abscense of antibiotic. Specifically, it counts, from xlsx
files obtained by parsing the breseq output

"""

import pandas as pd
import os
import pprint as pp

NJDirName = "/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_NJ_by_group" # Directory with NJ xlsx files for each condition
MutDirName = "/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_mutations_by_group" # Directory with mutations xlsx files for each condition

NJfiles = os.listdir(NJDirName)
Mutfiles = os.listdir(MutDirName)

Dic_Mutfiles = {}
strain_list = []
Dic_NJfiles = {}

for Mutfile in Mutfiles:
    Mutf_whole = MutDirName+"/"+Mutfile
    strain = Mutfile.split("_")[0]
    if strain not in strain_list:
        strain_list.append(strain)
    if "_nopOXA_mutation" in Mutf_whole: # just no pOXA condition
        #if "CF12" in Mutf_whole: #just CF12 to make test
        Mut_nopOXA = pd.read_excel(Mutf_whole)
        # print(Mut_nopOXA)
        annot_list = list(Mut_nopOXA["Annotation"])
        contig_list = list(Mut_nopOXA["Seq.ID"])
        mut_list = list(Mut_nopOXA["Mutation"])
        # print(annot_list)
        cols = [12, 13, 14]
        for col in cols:
            interg_nopOXA = 0
            nonsyn_nopOXA = 0
            syn_nopOXA = 0
            snps_nopOXA = 0
            indels_nopOXA = 0
            for i in range(len(Mut_nopOXA)):
                if list(Mut_nopOXA.iloc[:,col])[i] != 0:
                    if contig_list[i] == 1 or contig_list[i] == "contig_1": # if the mutation is in the chromosome
                        p_annot = annot_list[i].split("(")
                        # print(p_annot)
                        # print(p_annot[0][0], p_annot[0][-1])
                        if p_annot[0] == "intergenic ":
                            interg_nopOXA += 1
                        if p_annot[0] == "coding ":
                            nonsyn_nopOXA +=1
                        if p_annot[0] != "coding " and p_annot[0] != "intergenic " and p_annot[0][0] != p_annot[0][-1]:
                            nonsyn_nopOXA += 1
                        if p_annot[0][0] == p_annot[0][-1]:
                            syn_nopOXA += 1
                        if "+" in mut_list[i] or "Δ" in mut_list[i]:
                            indels_nopOXA += 1
                        if "→" in mut_list[i]:
                            snps_nopOXA += 1
        
            # dictionary with strain+condition as key and intergenic, nonsynonimous and synonimous mutations as values        
            Dic_Mutfiles[strain+"_nopOXA_"+str(col)] = [interg_nopOXA, nonsyn_nopOXA, syn_nopOXA, indels_nopOXA, snps_nopOXA]
            # print(Dic_Mutfiles)
        
    if "_pOXA_Ab_mutation" in Mutf_whole: # pOXA and Ab condition
        #if "CF12" in Mutf_whole: #just CF12 to make test
        Mut_pOXA_Ab = pd.read_excel(Mutf_whole)
        # print(Mut_pOXA)
        strain = Mutfile.split("_")[0]
        strain_list.append(strain)
        annot_list = list(Mut_pOXA_Ab["Annotation"])
        contig_list = list(Mut_pOXA_Ab["Seq.ID"])
        mut_list = list(Mut_pOXA_Ab["Mutation"])
        # print(annot_list)
        cols = [15, 16, 17]
        for col in cols:
            interg_pOXA_Ab = 0
            nonsyn_pOXA_Ab = 0
            syn_pOXA_Ab = 0
            snps_pOXA_Ab = 0
            indels_pOXA_Ab = 0
            for i in range(len(Mut_pOXA_Ab)):
                if list(Mut_pOXA_Ab.iloc[:,col])[i] != 0 and list(Mut_pOXA_Ab.iloc[:,8])[i] == 0 and list(Mut_pOXA_Ab.iloc[:,9])[i] == 0 and list(Mut_pOXA_Ab.iloc[:,10])[i] == 0:
                    if contig_list[i] == 1 or contig_list[i] == "contig_1": # if the mutation is in the chromosome
                        p_annot = str(annot_list[i]).split("(")
                        # print(p_annot)
                        # print(p_annot[0][0], p_annot[0][-1])
                        if p_annot[0] == "intergenic ":
                            interg_pOXA_Ab += 1
                        if p_annot[0] == "coding ":
                            nonsyn_pOXA_Ab +=1
                        if p_annot[0] != "coding " and p_annot[0] != "intergenic " and p_annot[0][0] != p_annot[0][-1]:
                            nonsyn_pOXA_Ab += 1
                        if p_annot[0][0] == p_annot[0][-1]:
                            syn_pOXA_Ab += 1
                        if "+" in mut_list[i] or "Δ" in mut_list[i]:
                            indels_pOXA_Ab += 1
                        if "→" in mut_list[i]:
                            snps_pOXA_Ab += 1
                
            # dictionary with strain+condition as key and intergenic, nonsynonimous and synonimous mutations as values
            Dic_Mutfiles[strain+"_pOXA_Ab_"+str(col)] = [interg_pOXA_Ab, nonsyn_pOXA_Ab, syn_pOXA_Ab, indels_pOXA_Ab, snps_pOXA_Ab]
            # print(Dic_Mutfiles)

# pp.pprint(Dic_Mutfiles)

# extract info from dictionary to obtain mutation spectrum of the chromosome in each strain
strain_list = list(dict.fromkeys(strain_list))
st_cond_list = Dic_Mutfiles.keys()
#print(st_cond_list)

for st_cond in Dic_Mutfiles:
    #print(Dic_Mutfiles[st_cond])
    totalmuts = Dic_Mutfiles[st_cond][0] + Dic_Mutfiles [st_cond][1] + Dic_Mutfiles [st_cond][2]
    if Dic_Mutfiles[st_cond][2] != 0:
        ns_s_ratio = Dic_Mutfiles[st_cond][1]/Dic_Mutfiles[st_cond][2]
    else:
        ns_s_ratio = str(Dic_Mutfiles[st_cond][1]) + "/" + str(Dic_Mutfiles[st_cond][2])
    interg_perc = (Dic_Mutfiles[st_cond][0]/totalmuts) * 100
    if Dic_Mutfiles[st_cond][3] != 0:
        snps_indel_ratio = Dic_Mutfiles[st_cond][4]/Dic_Mutfiles[st_cond][3]
    else:
        snps_indel_ratio = str(Dic_Mutfiles[st_cond][4]) +"/" + str(Dic_Mutfiles[st_cond][3])
    Dic_Mutfiles[st_cond].append(totalmuts) # include total number of mutations to values of dictionary
    Dic_Mutfiles[st_cond].append(ns_s_ratio) # include nonsyn/syn ratio
    Dic_Mutfiles[st_cond].append(interg_perc) # include intergenic percentage of mutations
    Dic_Mutfiles[st_cond].append(snps_indel_ratio) # include snps/indel ratio

df_data = pd.DataFrame(Dic_Mutfiles)
# print(df_data)
column_names = ["CF12_nopOXA_12", "CF12_nopOXA_13", "CF12_nopOXA_14",
                "CF12_pOXA_Ab_15", "CF12_pOXA_Ab_16", "CF12_pOXA_Ab_17",
                "CF13_nopOXA_12", "CF13_nopOXA_13", "CF13_nopOXA_14",
                "CF13_pOXA_Ab_15", "CF13_pOXA_Ab_16", "CF13_pOXA_Ab_17",
                "H53_nopOXA_12", "H53_nopOXA_13", "H53_nopOXA_14",
                "H53_pOXA_Ab_15", "H53_pOXA_Ab_16", "H53_pOXA_Ab_17",
                "K091_nopOXA_12", "K091_nopOXA_13", "K091_nopOXA_14",
                "K091_pOXA_Ab_15", "K091_pOXA_Ab_16", "K091_pOXA_Ab_17",
                "K147_nopOXA_12", "K147_nopOXA_13", "K147_nopOXA_14",
                "K147_pOXA_Ab_15", "K147_pOXA_Ab_16", "K147_pOXA_Ab_17",
                "K153_nopOXA_12", "K153_nopOXA_13", "K153_nopOXA_14",
                "K153_pOXA_Ab_15", "K153_pOXA_Ab_16", "K153_pOXA_Ab_17",
                "K163_nopOXA_12", "K163_nopOXA_13", "K163_nopOXA_14",
                "K163_pOXA_Ab_15", "K163_pOXA_Ab_16", "K163_pOXA_Ab_17",
                "K209_nopOXA_12", "K209_nopOXA_13", "K209_nopOXA_14",
                "K209_pOXA_Ab_15", "K209_pOXA_Ab_16", "K209_pOXA_Ab_17"]

df_tmp = df_data.reindex(columns = column_names)
df_final = df_tmp.reindex([5, 6, 0, 1, 2, 7, 3, 4, 8])
df_final = df_final.rename(index = {0: "Intergenic", 1: "Nonsynonymous", 
                                    2: "Synonymous", 3: "Indels",
                                    4: "SNPs", 5: "Total mutations", 
                                    6: "Nonsynonymous/Synonymous", 
                                    7: "Percent intergenic mutations",
                                    8: "SNPs/Indels"})
df_final = df_final.transpose()
print(df_final)

#df_final.to_excel("mutation_spectrum_final.xlsx")