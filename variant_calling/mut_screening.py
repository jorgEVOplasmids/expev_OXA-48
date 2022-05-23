#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 12:12:09 2022

@author: jorge

This script takes a breseq excel output parsed with the breseq_index_parser.py script
and generates an excel file summarizing which genes are mutated in which condition and
with what frequency (freq of mutations in each sample is summed for the same gene to avoid
duplications in the table). This is done for the strain indicated in the strain variable.

"""

import os
import pandas as pd
import seaborn as sb
from matplotlib import pyplot

DirName = "/home/jorge/Documents/TFM/clinical_samples/R_analysis/xlsx_breseqs"

files = os.listdir(DirName)

list_count_genes = []

strain = "K209"

for file in files:
    file_whole = DirName + "/" + file
    if strain in file: # for test
        xlsx = pd.ExcelFile(file_whole)
        muts = pd.read_excel(xlsx, "mutations")
        NJ = pd.read_excel(xlsx, "JC full")
        # For the mutations sheet, we'll collect lists of mutated genes without and with Ab
        # First of all, we define the necessary variables for this
        genes = muts["Gene"]
        contig = muts["Seq ID"]
        mut_position = muts["Position"]
        mut_pos_to_drop = []
        sample_pOXA_1 = muts.iloc[:,8]
        sample_pOXA_2 = muts.iloc[:,9]
        sample_pOXA_3 = muts.iloc[:,10]
        sample_pOXA_Ab_1 = muts.iloc[:,15]
        sample_pOXA_Ab_2 = muts.iloc[:,16]
        sample_pOXA_Ab_3 = muts.iloc[:,17]
        sample_pOXA_Anc = muts.iloc[:,7]
        sample_nopOXA_Anc = muts.iloc[:,11]
        sample_nopOXA_1 = muts.iloc[:,12]
        sample_nopOXA_2 = muts.iloc[:,13]
        sample_nopOXA_3 = muts.iloc[:,14]
        genes_pOXA_1 = {}
        genes_pOXA_2 = {}
        genes_pOXA_3 = {}
        genes_pOXA_Ab_1 = {}
        genes_pOXA_Ab_2 = {}
        genes_pOXA_Ab_3 = {}
        genes_nopOXA_1 = {}
        genes_nopOXA_2 = {}
        genes_nopOXA_3 = {}
        # Firstly, get a list of the targeted positions in the chromosome in the ancestors
        for i in range(len(muts)):
            if sample_pOXA_Anc[i] != 0 or sample_nopOXA_Anc[i] != 0:
                mut_pos_to_drop.append(mut_position[i])
        # Then we parse the xlsx and get the genes mutated in each sample
        for i in range(len(muts)):
            if sample_pOXA_1[i] != 0 and sample_pOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_pOXA_1:
                    genes_pOXA_1[genes[i]] = sample_pOXA_1[i]
                else:
                    genes_pOXA_1[genes[i]] += sample_pOXA_1[i]
            if sample_pOXA_2[i] != 0 and sample_pOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_pOXA_2:
                    genes_pOXA_2[genes[i]] = sample_pOXA_2[i]
                else:
                    genes_pOXA_2[genes[i]] += sample_pOXA_2[i]
            if sample_pOXA_3[i] != 0 and sample_pOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_pOXA_3:
                    genes_pOXA_3[genes[i]] = sample_pOXA_3[i]
                else:
                    genes_pOXA_3[genes[i]] += sample_pOXA_3[i]
            if sample_pOXA_Ab_1[i] != 0 and sample_pOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_pOXA_Ab_1:
                    genes_pOXA_Ab_1[genes[i]] = sample_pOXA_Ab_1[i]
                else:
                    genes_pOXA_Ab_1[genes[i]] += sample_pOXA_Ab_1[i]
            if sample_pOXA_Ab_2[i] != 0 and sample_pOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_pOXA_Ab_2:
                    genes_pOXA_Ab_2[genes[i]] = sample_pOXA_Ab_2[i]
                else:
                    genes_pOXA_Ab_2[genes[i]] += sample_pOXA_Ab_2[i]
            if sample_pOXA_Ab_3[i] != 0 and sample_pOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_pOXA_Ab_3:
                    genes_pOXA_Ab_3[genes[i]] = sample_pOXA_Ab_3[i]
                else:
                    genes_pOXA_Ab_3[genes[i]] += sample_pOXA_Ab_3[i]
            if sample_nopOXA_1[i] != 0 and sample_nopOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_nopOXA_1:
                    genes_nopOXA_1[genes[i]] = sample_nopOXA_1[i]
                else:
                    genes_nopOXA_1[genes[i]] += sample_nopOXA_1[i]
            if sample_nopOXA_2[i] != 0 and sample_nopOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_nopOXA_2:
                    genes_nopOXA_2[genes[i]] = sample_nopOXA_2[i]
                else:
                    genes_nopOXA_2[genes[i]] += sample_nopOXA_2[i]
            if sample_nopOXA_3[i] != 0 and sample_nopOXA_Anc[i] == 0 and mut_position[i] not in mut_pos_to_drop:
                if genes[i] not in genes_nopOXA_3:
                    genes_nopOXA_3[genes[i]] = sample_nopOXA_3[i]
                else:
                    genes_nopOXA_3[genes[i]] += sample_nopOXA_3[i]
        #list_count_genes.append(dict(collections.Counter(genes_pOXA_1)))
        #list_count_genes.append(dict(collections.Counter(genes_pOXA_2)))
        #list_count_genes.append(dict(collections.Counter(genes_pOXA_3)))
        #list_count_genes.append(dict(collections.Counter(genes_pOXA_Ab_1)))
        #list_count_genes.append(dict(collections.Counter(genes_pOXA_Ab_2)))
        #list_count_genes.append(dict(collections.Counter(genes_pOXA_Ab_3)))
        
        # For the NJ sheet we have to do some extra tasks
        # Define the necessary variables
        pos_to_drop = []
        NJ_contig = NJ["Seq ID 1"]
        NJ_position = NJ["Position 1"]
        NJ_gene1 = NJ["Gene 1"]
        NJ_gene2 = NJ["Gene 2"]
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
        # Firstly, get a list of the targeted positions in the chromosome in the ancestors
        for i in range(len(NJ)):
            if NJ_pOXA_Anc[i] != 0 or NJ_nopOXA_Anc[i] != 0:
                pos_to_drop.append(NJ_position[i])
        # Now, get the genes targeted by the new junctions
        # If the gene is targeted by both ends of the NJ it is only counted once
        for i in range(len(NJ)):
            if NJ_pOXA_1[i] != 0 and NJ_pOXA_1[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_pOXA_1:
                        genes_pOXA_1[NJ_gene1[i]] = NJ_pOXA_1[i]
                    else:
                        genes_pOXA_1[NJ_gene1[i]] += NJ_pOXA_1[i]
                    if NJ_gene2[i] not in genes_pOXA_1:
                        genes_pOXA_1[NJ_gene2[i]] = NJ_pOXA_1[i]
                    else:
                        genes_pOXA_1[NJ_gene2[i]] += NJ_pOXA_1[i]
                else:
                    if NJ_gene1[i] not in genes_pOXA_1:
                        genes_pOXA_1[NJ_gene1[i]] = NJ_pOXA_1[i]
                    else:
                        genes_pOXA_1[NJ_gene1[i]] += NJ_pOXA_1[i]
            if NJ_pOXA_2[i] != 0 and NJ_pOXA_2[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_pOXA_2:
                        genes_pOXA_2[NJ_gene1[i]] = NJ_pOXA_2[i]
                    else:
                        genes_pOXA_2[NJ_gene1[i]] += NJ_pOXA_2[i]
                    if NJ_gene2[i] not in genes_pOXA_2:
                        genes_pOXA_2[NJ_gene2[i]] = NJ_pOXA_2[i]
                    else:
                        genes_pOXA_2[NJ_gene2[i]] += NJ_pOXA_2[i]
                else:
                    if NJ_gene1[i] not in genes_pOXA_2:
                        genes_pOXA_2[NJ_gene1[i]] = NJ_pOXA_2[i]
                    else:
                        genes_pOXA_2[NJ_gene1[i]] += NJ_pOXA_2[i]
            if NJ_pOXA_3[i] != 0 and NJ_pOXA_3[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_pOXA_3:
                        genes_pOXA_3[NJ_gene1[i]] = NJ_pOXA_3[i]
                    else:
                        genes_pOXA_3[NJ_gene1[i]] += NJ_pOXA_3[i]
                    if NJ_gene2[i] not in genes_pOXA_3:
                        genes_pOXA_3[NJ_gene2[i]] = NJ_pOXA_3[i]
                    else:
                        genes_pOXA_3[NJ_gene2[i]] += NJ_pOXA_3[i]
                else:
                    if NJ_gene1[i] not in genes_pOXA_3:
                        genes_pOXA_3[NJ_gene1[i]] = NJ_pOXA_3[i]
                    else:
                        genes_pOXA_3[NJ_gene1[i]] += NJ_pOXA_3[i]
            if NJ_pOXA_Ab_1[i] != 0 and NJ_pOXA_Ab_1[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_pOXA_Ab_1:
                        genes_pOXA_Ab_1[NJ_gene1[i]] = NJ_pOXA_Ab_1[i]
                    else:
                        genes_pOXA_Ab_1[NJ_gene1[i]] += NJ_pOXA_Ab_1[i]
                    if NJ_gene2[i] not in genes_pOXA_Ab_1:
                        genes_pOXA_Ab_1[NJ_gene2[i]] = NJ_pOXA_Ab_1[i]
                    else:
                        genes_pOXA_Ab_1[NJ_gene2[i]] += NJ_pOXA_Ab_1[i]
                else:
                    if NJ_gene1[i] not in genes_pOXA_Ab_1:
                        genes_pOXA_Ab_1[NJ_gene1[i]] = NJ_pOXA_Ab_1[i]
                    else:
                        genes_pOXA_Ab_1[NJ_gene1[i]] += NJ_pOXA_Ab_1[i]
            if NJ_pOXA_Ab_2[i] != 0 and NJ_pOXA_Ab_2[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_pOXA_Ab_2:
                        genes_pOXA_Ab_2[NJ_gene1[i]] = NJ_pOXA_Ab_2[i]
                    else:
                        genes_pOXA_Ab_2[NJ_gene1[i]] += NJ_pOXA_Ab_2[i]
                    if NJ_gene2[i] not in genes_pOXA_Ab_2:
                        genes_pOXA_Ab_2[NJ_gene2[i]] = NJ_pOXA_Ab_2[i]
                    else:
                        genes_pOXA_Ab_2[NJ_gene2[i]] += NJ_pOXA_Ab_2[i]
                else:
                    if NJ_gene1[i] not in genes_pOXA_Ab_2:
                        genes_pOXA_Ab_2[NJ_gene1[i]] = NJ_pOXA_Ab_2[i]
                    else:
                        genes_pOXA_Ab_2[NJ_gene1[i]] += NJ_pOXA_Ab_2[i]
            if NJ_pOXA_Ab_3[i] != 0 and NJ_pOXA_Ab_3[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_pOXA_Ab_3:
                        genes_pOXA_Ab_3[NJ_gene1[i]] = NJ_pOXA_Ab_3[i]
                    else:
                        genes_pOXA_Ab_3[NJ_gene1[i]] += NJ_pOXA_Ab_3[i]
                    if NJ_gene2[i] not in genes_pOXA_Ab_3:
                        genes_pOXA_Ab_3[NJ_gene2[i]] = NJ_pOXA_Ab_3[i]
                    else:
                        genes_pOXA_Ab_3[NJ_gene2[i]] += NJ_pOXA_Ab_3[i]
                else:
                    if NJ_gene1[i] not in genes_pOXA_Ab_3:
                        genes_pOXA_Ab_3[NJ_gene1[i]] = NJ_pOXA_Ab_3[i]
                    else:
                        genes_pOXA_Ab_3[NJ_gene1[i]] += NJ_pOXA_Ab_3[i]
            if NJ_nopOXA_1[i] != 0 and NJ_nopOXA_1[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_nopOXA_1:
                        genes_nopOXA_1[NJ_gene1[i]] = NJ_nopOXA_1[i]
                    else:
                        genes_nopOXA_1[NJ_gene1[i]] += NJ_nopOXA_1[i]
                    if NJ_gene2[i] not in genes_nopOXA_1:
                        genes_nopOXA_1[NJ_gene2[i]] = NJ_nopOXA_1[i]
                    else:
                        genes_nopOXA_1[NJ_gene2[i]] += NJ_nopOXA_1[i]
                else:
                    if NJ_gene1[i] not in genes_nopOXA_1:
                        genes_nopOXA_1[NJ_gene1[i]] = NJ_nopOXA_1[i]
                    else:
                        genes_nopOXA_1[NJ_gene1[i]] += NJ_nopOXA_1[i]
            if NJ_nopOXA_2[i] != 0 and NJ_nopOXA_2[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_nopOXA_2:
                        genes_nopOXA_2[NJ_gene1[i]] = NJ_nopOXA_2[i]
                    else:
                        genes_nopOXA_2[NJ_gene1[i]] += NJ_nopOXA_2[i]
                    if NJ_gene2[i] not in genes_nopOXA_2:
                        genes_nopOXA_2[NJ_gene2[i]] = NJ_nopOXA_2[i]
                    else:
                        genes_nopOXA_2[NJ_gene2[i]] += NJ_nopOXA_2[i]
                else:
                    if NJ_gene1[i] not in genes_nopOXA_2:
                        genes_nopOXA_2[NJ_gene1[i]] = NJ_nopOXA_2[i]
                    else:
                        genes_nopOXA_2[NJ_gene1[i]] += NJ_nopOXA_2[i]
            if NJ_nopOXA_3[i] != 0 and NJ_nopOXA_3[i] != "nan" and NJ_position[i] not in pos_to_drop:
                if NJ_gene1[i] != NJ_gene2[i]:
                    if NJ_gene1[i] not in genes_nopOXA_3:
                        genes_nopOXA_3[NJ_gene1[i]] = NJ_nopOXA_3[i]
                    else:
                        genes_nopOXA_3[NJ_gene1[i]] += NJ_nopOXA_3[i]
                    if NJ_gene2[i] not in genes_nopOXA_3:
                        genes_nopOXA_3[NJ_gene2[i]] = NJ_nopOXA_3[i]
                    else:
                        genes_nopOXA_3[NJ_gene2[i]] += NJ_nopOXA_3[i]
                else:
                    if NJ_gene1[i] not in genes_nopOXA_3:
                        genes_nopOXA_3[NJ_gene1[i]] = NJ_nopOXA_3[i]
                    else:
                        genes_nopOXA_3[NJ_gene1[i]] += NJ_nopOXA_3[i]
                        
        #genes_pOXA_Ab_1 = {}
        #genes_pOXA_3 = {}
        list_count_genes.append(genes_pOXA_1)
        list_count_genes.append(genes_pOXA_2)
        list_count_genes.append(genes_pOXA_3)
        list_count_genes.append(genes_pOXA_Ab_1)
        list_count_genes.append(genes_pOXA_Ab_2)
        list_count_genes.append(genes_pOXA_Ab_3)
        list_count_genes.append(genes_nopOXA_1)
        list_count_genes.append(genes_nopOXA_2)
        list_count_genes.append(genes_nopOXA_3)
        
        
# print(list_count_genes)
muts = pd.DataFrame(list_count_genes)
muts = muts.T
muts = muts.fillna(0)
muts = muts.rename(columns = ({0:"pOXA-48 1", 1:"pOXA-48 2", 2: "pOXA-48 3",
                               3: "pOXA-48+Ab 1", 4: "pOXA-48+Ab 2", 5: "pOXA-48+Ab 3",
                               6: "no pOXA-48 1", 7: "no pOXA-48 2", 8: "no pOXA-48 3"}))

# muts = muts.sort_index(axis=0, ascending=True)
#print(muts.shape)
#print(muts)

# Create heatmap of mutations by gene and by replicate

pyplot.figure(figsize=(9, 34)) # dimensions of plot may vary depending on the number of genes mutated
heatmap = sb.heatmap(muts, cmap="BuPu", xticklabels=True, yticklabels=True, vmax = 100)
heatmap.set_xticklabels(heatmap.get_xticklabels(),rotation = 60, horizontalalignment='right')
heatmap.set_title("Heatmap of "+str(strain)+" Mutated Genes", fontsize = 22)


# Save the table with all the mutations as excel

muts.to_excel("/home/jorge/Documents/TFM/clinical_samples/R_analysis/all_mutations_heatmaps/"+strain+"_all_muts.xlsx")
