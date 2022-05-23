# Code for plotting plasmid copy number of clinical strains samples calculated from the coverage 

setwd("~/Documents/TFM/clinical_samples/R_analysis")

library(ggplot2)
library(ggrepel)
library(lme4)
library(nlme)
library(geoR)
library(dplyr)
library(afex)
library(ggpubr)
library(Rcmdr)
library(multcomp)

# open tables with plasmid copy number information

pcn_C232_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/C232-pcn.tsv", header=TRUE)
pcn_CF12_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/CF12-pcn.tsv", header=TRUE)
pcn_CF13_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/CF13-pcn.tsv", header=TRUE)
pcn_H53_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/H53-pcn.tsv", header=TRUE)
pcn_K091_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/K091-pcn.tsv", header=TRUE)
pcn_K147_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/K147-pcn.tsv", header=TRUE)
pcn_K153_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/K153-pcn.tsv", header=TRUE)
pcn_K163_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/K163-pcn.tsv", header=TRUE)
pcn_K209_table <- read.table(file="~/Documents/TFM/clinical_samples/plasmid_copy_number/K209-pcn.tsv", header=TRUE)

# open files with plasmid information from each sample

C232_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/C232/C232-contigs-reformat.txt")
CF12_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/CF12/CF12-contigs-reformat.txt")
CF13_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/CF13/CF13-contigs-reformat.txt")
H53_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/H53/H53-contigs-reformat.txt")
K091_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/K091/K091-contigs-reformat.txt")
K147_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/K147/K147-contigs-reformat.txt")
K153_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/K153/K153-contigs-reformat.txt")
K163_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/K163/K163-contigs-reformat.txt")
K209_plasmids <- read.table(file="~/Documents/TFM/clinical_samples/anvio/K209/K209-contigs-reformat.txt")

# assign pOXA-48 to correspondant contig (of length of 65 Kb approx)

# C232_plasmids
pcn_C232_table$p_contig[pcn_C232_table$p_contig=="c_000000000003"] <- "pOXA-48"
# pcn_C232_table
# CF12_plasmids
pcn_CF12_table$p_contig[pcn_CF12_table$p_contig=="c_000000000002"] <- "pOXA-48"
# CF13_plasmids
pcn_CF13_table$p_contig[pcn_CF13_table$p_contig=="c_000000000002"] <- "pOXA-48"
# H53_plasmids
pcn_H53_table$p_contig[pcn_H53_table$p_contig=="c_000000000003"] <- "pOXA-48"
# K091_plasmids
pcn_K091_table$p_contig[pcn_K091_table$p_contig=="c_000000000003"] <- "pOXA-48"
# K147_plasmids
pcn_K147_table$p_contig[pcn_K147_table$p_contig=="c_000000000003"] <- "pOXA-48"
# K153_plasmids
pcn_K153_table$p_contig[pcn_K153_table$p_contig=="c_000000000002"] <- "pOXA-48"
# K163_plasmids
pcn_K163_table$p_contig[pcn_K163_table$p_contig=="c_000000000003"] <- "pOXA-48"
# K209_plasmids
pcn_K209_table$p_contig[pcn_K209_table$p_contig=="c_000000000004"] <- "pOXA-48"

# and the rest of plasmids

pcn_C232_table$p_contig[pcn_C232_table$p_contig=="c_000000000002"] <- "IncY_1"
pcn_C232_table$p_contig[pcn_C232_table$p_contig=="c_000000000004"] <- "pTORI_length=47289"
pcn_C232_table$p_contig[pcn_C232_table$p_contig=="c_000000000005"] <- "Col440I_1"
pcn_C232_table$p_contig[pcn_C232_table$p_contig=="c_000000000006"] <- "Col156_1"
pcn_C232_table$p_contig[pcn_C232_table$p_contig=="c_000000000007"] <- "Col156_2"
pcn_CF12_table$p_contig[pcn_CF12_table$p_contig=="c_000000000003"] <- "IncX3_1"
pcn_CF12_table$p_contig[pcn_CF12_table$p_contig=="c_000000000004"] <- "IncR_1"
pcn_CF12_table$p_contig[pcn_CF12_table$p_contig=="c_000000000005"] <- "ColE10_1"
pcn_CF12_table$p_contig[pcn_CF12_table$p_contig=="c_000000000006"] <- "ColRNAI_1"
pcn_CF13_table$p_contig[pcn_CF13_table$p_contig=="c_000000000003"] <- "IncR_1"
pcn_CF13_table$p_contig[pcn_CF13_table$p_contig=="c_000000000004"] <- "ColRNAI_1"
pcn_H53_table$p_contig[pcn_H53_table$p_contig=="c_000000000002"] <- "IncFIB(K)_1_Kpn3"
pcn_K091_table$p_contig[pcn_K091_table$p_contig=="c_000000000002"] <- "IncFIA(HI1)_1_HI1/IncFII_1_pKP91"
pcn_K147_table$p_contig[pcn_K147_table$p_contig=="c_000000000002"] <- "IncFIB(K)_1_Kpn3/IncFII_1_pKP91"
pcn_K153_table$p_contig[pcn_K153_table$p_contig=="c_000000000003"] <- "IncR_1"
pcn_K163_table$p_contig[pcn_K163_table$p_contig=="c_000000000002"] <- "IncFIB(K)_1_Kpn3/IncFII_1_pKP91"
pcn_K163_table$p_contig[pcn_K163_table$p_contig=="c_000000000004"] <- "ColRNAI_1"
pcn_K209_table$p_contig[pcn_K209_table$p_contig=="c_000000000002"] <- "IncHI1B_1_pNDM-MAR"
pcn_K209_table$p_contig[pcn_K209_table$p_contig=="c_000000000003"] <- "IncFIA(HI1)_1_HI1/IncFII_1_pKP91"

# assign group values to each sample in each strain (p, p+Ab, *p)

C232_OXA <- pcn_C232_table[pcn_C232_table$p_contig=="pOXA-48",]
C232_OXA$group <- c(rep("pOXA+Ab",3), rep("pOXA", 4))
C232_OXA$group <- factor(C232_OXA$group, levels = c("*pOXA" ,"pOXA", "pOXA+Ab"))
plot_C232_OXA <- C232_OXA[C232_OXA$group!="*pOXA",]
CF12_OXA <- pcn_CF12_table[pcn_CF12_table$p_contig=="pOXA-48",]
CF12_OXA$group <- c(rep("pOXA", 4), rep("*pOXA", 4), rep("pOXA+Ab", 3))
CF12_OXA$group <- factor(CF12_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_CF12_OXA <- CF12_OXA[CF12_OXA$group!="*pOXA",]
CF13_OXA <- pcn_CF13_table[pcn_CF13_table$p_contig=="pOXA-48",]
CF13_OXA$group <- c(rep("pOXA", 4), rep("*pOXA", 3), rep("pOXA+Ab", 3))
CF13_OXA$group <- factor(CF13_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_CF13_OXA <- CF13_OXA[CF13_OXA$group!="*pOXA",]
H53_OXA <- pcn_H53_table[pcn_H53_table$p_contig=="pOXA-48",]
H53_OXA$group <- c(rep("pOXA", 4), rep("*pOXA", 4), rep("pOXA+Ab", 3))
H53_OXA$group <- factor(H53_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_H53_OXA <- H53_OXA[H53_OXA$group!="*pOXA",]
K091_OXA <- pcn_K091_table[pcn_K091_table$p_contig=="pOXA-48",]
K091_OXA$group <- c(rep("pOXA", 4), rep("*pOXA", 4), rep("pOXA+Ab", 3))
K091_OXA$group <- factor(K091_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_K091_OXA <- K091_OXA[K091_OXA$group!="*pOXA",]
K147_OXA <- pcn_K147_table[pcn_K147_table$p_contig=="pOXA-48",]
K147_OXA$group <- c(rep("pOXA", 4), rep("*pOXA", 4), rep("pOXA+Ab", 3))
K147_OXA$group <- factor(K147_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_K147_OXA <- K147_OXA[K147_OXA$group!="*pOXA",]
K153_OXA <- pcn_K153_table[pcn_K153_table$p_contig=="pOXA-48",]
K153_OXA$group <- c(rep("pOXA+Ab", 3), rep("pOXA", 4), rep("*pOXA", 4))
K153_OXA$group <- factor(K153_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_K153_OXA <- K153_OXA[K153_OXA$group!="*pOXA",]
K163_OXA <- pcn_K163_table[pcn_K163_table$p_contig=="pOXA-48",]
K163_OXA$group <- c(rep("pOXA", 4), rep("*pOXA", 4), rep("pOXA+Ab", 3))
K163_OXA$group <- factor(K163_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_K163_OXA <- K163_OXA[K163_OXA$group!="*pOXA",]
K209_OXA <- pcn_K209_table[pcn_K209_table$p_contig=="pOXA-48",]
K209_OXA$group <- c(rep("pOXA", 4), rep(c("pOXA+Ab", "*pOXA"), 3), "*pOXA")
K209_OXA$group <- factor(K209_OXA$group, levels = c("*pOXA", "pOXA", "pOXA+Ab"))
plot_K209_OXA <- K209_OXA[K209_OXA$group!="*pOXA",]

# for all the plasmids in each sample, we have to shape the pcn_$strain_table and add the sample group (p, p+Ab, *p, p-Anc, *p-Anc)

pcn_C232_table$group <- c(rep("pOXA+Ab",18), rep("pOXA", 24))
pcn_CF12_table$group <- c(rep("pOXA", 15), rep("pOXA-Anc", 5), rep("*pOXA", 15), rep("*pOXA-Anc", 5), rep("pOXA+Ab", 15))
pcn_CF13_table$group <- c(rep("pOXA", 9), rep("pOXA-Anc", 3), rep("*pOXA", 6), rep("*pOXA-Anc", 3), rep("pOXA+Ab", 9))
pcn_H53_table$group <- c(rep("pOXA", 6), rep("pOXA-Anc", 2), rep("*pOXA", 6), rep("*pOXA-Anc", 2), rep("pOXA+Ab", 6))
pcn_K091_table$group <- c(rep("pOXA", 6), rep("pOXA-Anc", 2), rep("*pOXA", 6), rep("*pOXA-Anc", 2), rep("pOXA+Ab", 6))
pcn_K147_table$group <- c(rep("pOXA", 6), rep("pOXA-Anc", 2), rep("*pOXA", 6), rep("*pOXA-Anc", 2), rep("pOXA+Ab", 6))
pcn_K153_table$group <- c(rep("pOXA+Ab", 6), rep("pOXA", 6), rep("pOXA-Anc", 2), rep("*pOXA", 6), rep("*pOXA-Anc", 2))
pcn_K163_table$group <- c(rep("pOXA", 9), rep("pOXA-Anc", 3), rep("*pOXA", 9), rep("*pOXA-Anc", 3), rep("pOXA+Ab", 9))
pcn_K209_table$group <- c(rep("pOXA", 9), rep("pOXA-Anc", 3), rep(c(rep("pOXA+Ab", 3), rep("*pOXA", 3)), 3), rep("*pOXA-Anc", 3))

# reordering experimental groups for plots

pcn_CF12_table$group <- factor(pcn_CF12_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_CF13_table$group <- factor(pcn_CF13_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_H53_table$group <- factor(pcn_H53_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_K091_table$group <- factor(pcn_K091_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_K147_table$group <- factor(pcn_K147_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_K153_table$group <- factor(pcn_K153_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_K163_table$group <- factor(pcn_K163_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))
pcn_K209_table$group <- factor(pcn_K209_table$group, levels = c("*pOXA-Anc", "*pOXA", "pOXA-Anc", "pOXA", "pOXA+Ab"))

# final plots for all plasmids

ggplot(pcn_C232_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of C232 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_C232_table$p_contig, scales = "free")
ggplot(pcn_CF12_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of CF12 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_CF12_table$p_contig, scales = "free")
ggplot(pcn_CF13_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of CF13 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_CF13_table$p_contig, scales = "free")
ggplot(pcn_H53_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of H53 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_H53_table$p_contig)
ggplot(pcn_K091_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of K091 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_K091_table$p_contig)
ggplot(pcn_K147_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of K147 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_K147_table$p_contig, scales = "free")
ggplot(pcn_K153_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of K153 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_K153_table$p_contig)
ggplot(pcn_K163_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of K163 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_K163_table$p_contig, scales = "free")
ggplot(pcn_K209_table, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) +  geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.3, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change") + ggtitle("Coverage F-change of K209 plasmids") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=22)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) + theme(axis.text.y = element_text(size = 12)) + facet_wrap(~pcn_K209_table$p_contig, scales = "free")

# Filter out other plasmids different to pOXA-48

sa_H53 <- pcn_H53_table %>%
  filter(p_contig == "pOXA-48", group != "*pOXA" & group != "*pOXA-Anc")

sa_K091 <- pcn_K091_table %>%
  filter(p_contig == "pOXA-48", group != "*pOXA" & group != "*pOXA-Anc")

sa_K147 <- pcn_K147_table %>%
  filter(p_contig == "pOXA-48", group != "*pOXA" & group != "*pOXA-Anc")

sa_K153 <- pcn_K153_table %>%
  filter(p_contig == "pOXA-48", group != "*pOXA" & group != "*pOXA-Anc")

sa_K163 <- pcn_K163_table %>%
  filter(p_contig == "pOXA-48", group != "*pOXA" & group != "*pOXA-Anc")

sa_K209 <- pcn_K209_table %>%
  filter(p_contig == "pOXA-48", group != "*pOXA" & group != "*pOXA-Anc")

# STATISTICAL ANALYSIS FOR ALL SAMPLES ANALISED AT THE SAME TIME (WITHOUT C232 and CF12/CF13)

# generation of data frame and visualization of boxplot

total_OXA <- Reduce(function(x, y) merge(x, y, all=TRUE), list(sa_H53, sa_K091, sa_K147, sa_K153, sa_K163, sa_K209))
ggplot(total_OXA, aes(group, median_ratio)) + geom_point(aes(x = group, y = median_ratio, col = group)) + geom_boxplot(aes(x = group, y = median_ratio, col = group), width=0.5, alpha = 0.2) + xlab("Sample group") + ylab("Sequencing depth F-change (median cov. plasmid/median cov. chr)") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12))

total_OXA_anc <- total_OXA %>%
  filter(group != "pOXA")

# Comparison between the three groups simoultaneously

anovapOXA <- aov(median_ratio ~ group, data = total_OXA)
summary(anovapOXA)
numSummary(total_OXA$median_ratio, groups = total_OXA$group, statistics = c("mean", "sd"))

Pairs <- glht(anovapOXA, linfct = mcp(group = "Tukey"))
summary(Pairs)
confint(Pairs)
cld(Pairs)
old.oma <- par(oma = c(0,7,0,0))
plot(confint(Pairs))
par(old.oma)

# build LME model and check likelihood, here comparing ancestors to evolved samples with the antibiotic
total_LME <- lmer(median_ratio ~ group + (1+group|Strain), data = total_OXA_anc, REML = FALSE)
summary(total_LME)
total_LME.null <- lmer(median_ratio ~ (1+group|Strain), data = total_OXA_anc, REML = FALSE)
summary(total_LME.null)
anova(total_LME, total_LME.null)

# as the comparison of the models is significant, we can now validate our model

# Stansardised residuals vs fitted values
a <- plot(total_LME, resid(., scaled = TRUE) ~ fitted(.), abline = 0, pch = 16, xlab="Fitted values",ylab="Standardised residuals", main = "Stansardised residuals vs fitted values")
# Normal Q-Q
b <- ggplot(data = total_OXA_anc, aes(sample = resid(total_LME))) + stat_qq() + stat_qq_line() + xlab("Theoretical quantiles") + ylab("Standardized residuals") + ggtitle("Normal Q-Q") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=18)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12))
shapiro.test(resid(total_LME))
ks.test(resid(total_LME), pnorm)
# Histogram for residuals distribution
c <- ggplot(data=total_OXA_anc, aes(x=resid(total_LME))) + 
  geom_histogram(aes(y=..density..), bins=20) + 
  xlab("Residuals") + ylab("Density") +
  ggtitle("Histogram of residuals distribution") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=18)) +
  stat_function(fun=dnorm, 
                args=list(mean=mean(resid(total_LME)), 
                          sd=sd(resid(total_LME))))
# Scale-Location plot
d <- plot(total_LME, sqrt(abs(resid(., scaled = TRUE))) ~ fitted(.), abline = 0, pch = 16, xlab="Fitted values",ylab="Standardised residuals", main = "Scale-location plot")

figure <- ggarrange(a, b, c,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure

# repeat the comparison but with comparing the ancestors to samples evolved without antibiotic

total_OXA_noab <- total_OXA %>%
  filter(group != "pOXA+Ab")

total_LME <- lmer(median_ratio ~ group + (1+group|Strain), data = total_OXA_noab, REML = FALSE)
summary(total_LME)
total_LME.null <- lmer(median_ratio ~ (1+group|Strain), data = total_OXA_noab, REML = FALSE)
summary(total_LME.null)
anova(total_LME, total_LME.null)

# as the comparison of the models is significant, we can now validate our model

# Stansardised residuals vs fitted values
a <- plot(total_LME, resid(., scaled = TRUE) ~ fitted(.), abline = 0, pch = 16, xlab="Fitted values",ylab="Standardised residuals", main = "Stansardised residuals vs fitted values")
# Normal Q-Q
b <- ggplot(data = total_OXA_noab, aes(sample = resid(total_LME))) + stat_qq() + stat_qq_line() + xlab("Theoretical quantiles") + ylab("Standardized residuals") + ggtitle("Normal Q-Q") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=18)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12))
shapiro.test(resid(total_LME))
ks.test(resid(total_LME), pnorm)
# Histogram for residuals distribution
c <- ggplot(data=total_OXA_noab, aes(x=resid(total_LME))) + 
  geom_histogram(aes(y=..density..), bins=20) + 
  xlab("Residuals") + ylab("Density") +
  ggtitle("Histogram of residuals distribution") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=18)) +
  stat_function(fun=dnorm, 
                args=list(mean=mean(resid(total_LME)), 
                          sd=sd(resid(total_LME))))
# Scale-Location plot
d <- plot(total_LME, sqrt(abs(resid(., scaled = TRUE))) ~ fitted(.), abline = 0, pch = 16, xlab="Fitted values",ylab="Standardised residuals", main = "Scale-location plot")

figure <- ggarrange(a, b, c,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure


# and now with no ancestors

total_OXA_noanc <- total_OXA %>%
  filter(group != "pOXA-Anc")

total_LME <- lmer(median_ratio ~ group + (1+group|Strain), data = total_OXA_noanc, REML = FALSE)
summary(total_LME)
total_LME.null <- lmer(median_ratio ~ (1+group|Strain), data = total_OXA_noanc, REML = FALSE)
summary(total_LME.null)
anova(total_LME, total_LME.null)

# as the comparison of the models is significant, we can now validate our model

# Stansardised residuals vs fitted values
a <- plot(total_LME, resid(., scaled = TRUE) ~ fitted(.), abline = 0, pch = 16, xlab="Fitted values",ylab="Standardised residuals", main = "Stansardised residuals vs fitted values")
# Normal Q-Q
b <- ggplot(data = total_OXA_noanc, aes(sample = resid(total_LME))) + stat_qq() + stat_qq_line() + xlab("Theoretical quantiles") + ylab("Standardized residuals") + ggtitle("Normal Q-Q") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=18)) + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12))
shapiro.test(resid(total_LME))
ks.test(resid(total_LME), pnorm)
# Histogram for residuals distribution
c <- ggplot(data=total_OXA_noanc, aes(x=resid(total_LME))) + 
  geom_histogram(aes(y=..density..), bins=20) + 
  xlab("Residuals") + ylab("Density") +
  ggtitle("Histogram of residuals distribution") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=18)) +
  stat_function(fun=dnorm, 
                args=list(mean=mean(resid(total_LME)), 
                          sd=sd(resid(total_LME))))
# Scale-Location plot
d <- plot(total_LME, sqrt(abs(resid(., scaled = TRUE))) ~ fitted(.), abline = 0, pch = 16, xlab="Fitted values",ylab="Standardised residuals", main = "Scale-location plot")

figure <- ggarrange(a, b, c,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure
