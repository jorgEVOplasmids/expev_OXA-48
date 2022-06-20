setwd("~/Universidad/TFM/Curvas")

library(tidyverse)
library(readxl)
library(lattice)
library(deSolve)
library(growthrates)
library(dplyr)
library(gridExtra)
library(caTools)
library(tidyr)
library(flux)
library(ggrepel)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(car)
#library(multcomp)

## Statistical test for plasmid cost

# CF12
# pOXA-48 vs no pOXA-48

CF12_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = CF12_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF12_aov_data$AUC,CF12_aov_data$Plasmid)

# Compute anova

CF12.aov <- aov(CF12_aov_data$AUC ~ CF12_aov_data$Plasmid); CF12.aov
summary(CF12.aov)

# pOXA-48+Ab vs no pOXA-48

CF12_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

CF12_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

CF12_ab_aov_data <- rbind(CF12_aov_data_nopoxa, CF12_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = CF12_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF12_ab_aov_data$AUC,CF12_ab_aov_data$Plasmid)

# Compute anova

CF12.aov <- aov(CF12_ab_aov_data$AUC ~ CF12_ab_aov_data$Plasmid); CF12.aov
summary(CF12.aov)

# CF13
# pOXA-48 vs no pOXA-48

CF13_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = CF13_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF13_aov_data$AUC,CF13_aov_data$Plasmid)

# Compute anova

CF13.aov <- aov(CF13_aov_data$AUC ~ CF13_aov_data$Plasmid); CF13.aov
summary(CF13.aov)

# pOXA-48+Ab vs no pOXA-48

CF13_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

CF13_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

CF13_ab_aov_data <- rbind(CF13_aov_data_nopoxa, CF13_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = CF13_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF13_ab_aov_data$AUC,CF13_ab_aov_data$Plasmid)

# Compute anova

CF13.aov <- aov(CF13_ab_aov_data$AUC ~ CF13_ab_aov_data$Plasmid); CF13.aov
summary(CF13.aov)

# H53
# pOXA-48 vs no pOXA-48

H53_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = H53_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(H53_aov_data$AUC,H53_aov_data$Plasmid)

# Compute anova

H53.aov <- aov(H53_aov_data$AUC ~ H53_aov_data$Plasmid); H53.aov
summary(H53.aov)

# pOXA-48+Ab vs no pOXA-48

H53_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

H53_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

H53_ab_aov_data <- rbind(H53_aov_data_nopoxa, H53_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = H53_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(H53_ab_aov_data$AUC,H53_ab_aov_data$Plasmid)

# Compute anova

H53.aov <- aov(H53_ab_aov_data$AUC ~ H53_ab_aov_data$Plasmid); H53.aov
summary(H53.aov)

# K091
# pOXA-48 vs no pOXA-48

K091_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K091_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K091_aov_data$AUC,K091_aov_data$Plasmid)

# Compute anova

K091.aov <- aov(K091_aov_data$AUC ~ K091_aov_data$Plasmid); K091.aov
summary(K091.aov)

# pOXA-48+Ab vs no pOXA-48

K091_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

K091_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

K091_ab_aov_data <- rbind(K091_aov_data_nopoxa, K091_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K091_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K091_ab_aov_data$AUC,K091_ab_aov_data$Plasmid)

# Compute anova

K091.aov <- aov(K091_ab_aov_data$AUC ~ K091_ab_aov_data$Plasmid); K091.aov
summary(K091.aov)

# K147
# pOXA-48 vs no pOXA-48

K147_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K147_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K147_aov_data$AUC,K147_aov_data$Plasmid)

# Compute anova

K147.aov <- aov(K147_aov_data$AUC ~ K147_aov_data$Plasmid); K147.aov
summary(K147.aov)

# pOXA-48+Ab vs no pOXA-48

K147_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

K147_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

K147_ab_aov_data <- rbind(K147_aov_data_nopoxa, K147_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K147_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K147_ab_aov_data$AUC,K147_ab_aov_data$Plasmid)

# Compute anova

K147.aov <- aov(K147_ab_aov_data$AUC ~ K147_ab_aov_data$Plasmid); K147.aov
summary(K147.aov)

# K153
# pOXA-48 vs no pOXA-48

K153_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K153_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K153_aov_data$AUC,K153_aov_data$Plasmid)

# Compute anova

K153.aov <- aov(K153_aov_data$AUC ~ K153_aov_data$Plasmid); K153.aov
summary(K153.aov)

# pOXA-48+Ab vs no pOXA-48

K153_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

K153_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

K153_ab_aov_data <- rbind(K153_aov_data_nopoxa, K153_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K153_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K153_ab_aov_data$AUC,K153_ab_aov_data$Plasmid)

# Compute anova

K153.aov <- aov(K153_ab_aov_data$AUC ~ K153_ab_aov_data$Plasmid); K153.aov
summary(K153.aov)

# K163
# pOXA-48 vs no pOXA-48

K163_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K163_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K163_aov_data$AUC,K163_aov_data$Plasmid)

# Compute anova

K163.aov <- aov(K163_aov_data$AUC ~ K163_aov_data$Plasmid); K163.aov
summary(K163.aov)

# pOXA-48+Ab vs no pOXA-48

K163_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

K163_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

K163_ab_aov_data <- rbind(K163_aov_data_nopoxa, K163_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K163_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K163_ab_aov_data$AUC,K163_ab_aov_data$Plasmid)

# Compute anova

K163.aov <- aov(K163_ab_aov_data$AUC ~ K163_ab_aov_data$Plasmid); K163.aov
summary(K163.aov)

# K209
# pOXA-48 vs no pOXA-48

K209_aov_data <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & Day == "1" & Antibiotic == "NO")

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K209_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_aov_data$AUC,K209_aov_data$Plasmid)

# Compute anova

K209.aov <- aov(K209_aov_data$AUC ~ K209_aov_data$Plasmid); K209.aov
summary(K209.aov)

# pOXA-48+Ab vs no pOXA-48

K209_aov_data_nopoxa <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & Day == "1") %>%
  filter(Plasmid != "YES" & Antibiotic == "NO")

K209_aov_data_poxaab <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & Day == "1") %>%
  filter(Plasmid == "YES" & Antibiotic == "YES")

K209_ab_aov_data <- rbind(K209_aov_data_nopoxa, K209_aov_data_poxaab)

# Check normality assumption

model  <- lm(AUC ~ Plasmid, data = K209_ab_aov_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_ab_aov_data$AUC,K209_ab_aov_data$Plasmid)

# Compute anova

K209.aov <- aov(K209_ab_aov_data$AUC ~ K209_ab_aov_data$Plasmid); K209.aov
summary(K209.aov)

## Correct p values by FDR
# pOXA-48 vs no pOXA
plist_noAb <- c(0.746,0.968,0.0513,0.174,0.909,0.315,0.102,0.000117)
p.adjust(plist_noAb, method = "fdr", n = length(plist_noAb))

# pOXA-48+Ab vs no pOXA
plist_Ab <- c(0.00148,0.0675,0.00866,0.0295,0.648,0.228,0.252)
p.adjust(plist_Ab, method = "fdr", n = length(plist_Ab))

### Statistical analysis of fitness increase between days 1 and 15

## CF12
# no pOXA

CF12_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = CF12_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF12_fit_data$AUC,CF12_fit_data$Day)

# Compute anova

CF12.aov <- aov(CF12_fit_data$AUC ~ CF12_fit_data$Day); CF12.aov
summary(CF12.aov)

# pOXA-48

CF12_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = CF12_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF12_fit_data$AUC,CF12_fit_data$Day)

# Compute anova

CF12.aov <- aov(CF12_fit_data$AUC ~ CF12_fit_data$Day); CF12.aov
summary(CF12.aov)

# pOXA-48+Ab

CF12_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = CF12_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF12_fit_data$AUC,CF12_fit_data$Day)

# Compute anova

CF12.aov <- aov(CF12_fit_data$AUC ~ CF12_fit_data$Day); CF12.aov
summary(CF12.aov)

#----------------------------------------------------------------

## CF13
# no pOXA

CF13_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = CF13_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF13_fit_data$AUC,CF13_fit_data$Day)

# Compute anova

CF13.aov <- aov(CF13_fit_data$AUC ~ CF13_fit_data$Day); CF13.aov
summary(CF13.aov)

# pOXA-48

CF13_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = CF13_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF13_fit_data$AUC,CF13_fit_data$Day)

# Compute anova

CF13.aov <- aov(CF13_fit_data$AUC ~ CF13_fit_data$Day); CF13.aov
summary(CF13.aov)

# pOXA-48+Ab

CF13_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = CF13_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF13_fit_data$AUC,CF13_fit_data$Day)

# Compute anova

CF13.wilk <- wilcox.test(CF13_fit_data$AUC ~ CF13_fit_data$Day); CF13.wilk

CF13.aov <- aov(CF13_fit_data$AUC ~ CF13_fit_data$Day); CF13.aov
summary(CF13.aov)

# ----------------------------------------------------

## H53
# no pOXA

H53_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = H53_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(H53_fit_data$AUC,H53_fit_data$Day)

# Compute anova

H53.aov <- aov(H53_fit_data$AUC ~ H53_fit_data$Day); H53.aov
summary(H53.aov)

# pOXA-48

H53_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = H53_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(H53_fit_data$AUC,H53_fit_data$Day)

# Compute anova

H53.aov <- aov(H53_fit_data$AUC ~ H53_fit_data$Day); H53.aov
summary(H53.aov)

# pOXA-48+Ab

H53_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = H53_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(H53_fit_data$AUC,H53_fit_data$Day)

# Compute anova

H53.aov <- aov(H53_fit_data$AUC ~ H53_fit_data$Day); H53.aov
summary(H53.aov)

# ------------------------------------------------------------------

## K091
# no pOXA

K091_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K091_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K091_fit_data$AUC,K091_fit_data$Day)

# Compute anova

K091.aov <- aov(K091_fit_data$AUC ~ K091_fit_data$Day); K091.aov
summary(K091.aov)

# pOXA-48

K091_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K091_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K091_fit_data$AUC,K091_fit_data$Day)

# Compute anova

K091.aov <- aov(K091_fit_data$AUC ~ K091_fit_data$Day); K091.aov
summary(K091.aov)

# pOXA-48+Ab

K091_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K091_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K091_fit_data$AUC,K091_fit_data$Day)

# Compute anova

K091.aov <- aov(K091_fit_data$AUC ~ K091_fit_data$Day); K091.aov
summary(K091.aov)

# -----------------------------------------------------------------------------

## K147
# no pOXA

K147_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K147_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K147_fit_data$AUC,K147_fit_data$Day)

# Compute anova

K147.aov <- aov(K147_fit_data$AUC ~ K147_fit_data$Day); K147.aov
summary(K147.aov)

# pOXA-48

K147_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K147_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K147_fit_data$AUC,K147_fit_data$Day)

# Compute anova

K147.aov <- aov(K147_fit_data$AUC ~ K147_fit_data$Day); K147.aov
summary(K147.aov)

# pOXA-48+Ab

K147_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K147_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K147_fit_data$AUC,K147_fit_data$Day)

# Compute anova

K147.aov <- aov(K147_fit_data$AUC ~ K147_fit_data$Day); K147.aov
summary(K147.aov)

# --------------------------------------------------------------------

## K153
# no pOXA

K153_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K153_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K153_fit_data$AUC,K153_fit_data$Day)

# Compute anova

K153.aov <- aov(K153_fit_data$AUC ~ K153_fit_data$Day); K153.aov
summary(K153.aov)

# pOXA-48

K153_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K153_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K153_fit_data$AUC,K153_fit_data$Day)

# Compute anova

K153.aov <- aov(K153_fit_data$AUC ~ K153_fit_data$Day); K153.aov
summary(K153.aov)

# pOXA-48+Ab

K153_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K153_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K153_fit_data$AUC,K153_fit_data$Day)

# Compute anova

K153.aov <- aov(K153_fit_data$AUC ~ K153_fit_data$Day); K153.aov
summary(K153.aov)

# --------------------------------------------------------------

## K163
# no pOXA

K163_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K163_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K163_fit_data$AUC,K163_fit_data$Day)

# Compute anova

K163.aov <- aov(K163_fit_data$AUC ~ K163_fit_data$Day); K163.aov
summary(K163.aov)

# pOXA-48

K163_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K163_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K163_fit_data$AUC,K163_fit_data$Day)

# Compute anova

K163.wilk <- wilcox.test(K163_fit_data$AUC ~ K163_fit_data$Day); K163.wilk
summary(K163.wilk)

K163.aov <- aov(K163_fit_data$AUC ~ K163_fit_data$Day); K163.aov
summary(K163.aov)

# pOXA-48+Ab

K163_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K163_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K163_fit_data$AUC,K163_fit_data$Day)

# Compute anova

K163.wilk <- wilcox.test(K163_fit_data$AUC ~ K163_fit_data$Day); K163.wilk

K163.aov <- aov(K163_fit_data$AUC ~ K163_fit_data$Day); K163.aov
summary(K163.aov)

# ------------------------------------------------------------------------

## K209
# no pOXA

K209_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "NO")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K209_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_fit_data$AUC,K209_fit_data$Day)

# Compute anova

K209.aov <- aov(K209_fit_data$AUC ~ K209_fit_data$Day); K209.aov
summary(K209.aov)

# pOXA-48

K209_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "NO" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K209_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_fit_data$AUC,K209_fit_data$Day)

# Compute anova

K209.aov <- aov(K209_fit_data$AUC ~ K209_fit_data$Day); K209.aov
summary(K209.aov)

# pOXA-48+Ab

K209_fit_data <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & ( Day == "1" | Day == "15")) %>%
  filter(Antibiotic == "YES" & Plasmid == "YES")

# Check normality assumption

model  <- lm(AUC ~ Day, data = K209_fit_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_fit_data$AUC,K209_fit_data$Day)

# Compute anova

K209.aov <- aov(K209_fit_data$AUC ~ K209_fit_data$Day); K209.aov
summary(K209.aov)

# Correct p-values by FDR
# no pOXA
plist_noAb <- c(1.87e-05,2.17e-06,0.000409,1.52e-06,4.47e-05,0.000451,0.000564,0.0585)
p.adjust(plist_noAb, method = "fdr", n = length(plist_noAb))

# pOXA-48
plist_Ab <- c(0.00056,2.59e-05,3.3e-06,0.000274,4.09e-07,3.56e-05,0.002165,0.954)
p.adjust(plist_Ab, method = "fdr", n = length(plist_Ab))

# pOXA-48+Ab
plist_Ab <- c(9.75e-05,0.002165,0.000266,1.43e-05,2.99e-05,1.78e-05,0.002165,0.3)
p.adjust(plist_Ab, method = "fdr", n = length(plist_Ab))

### Compare fitness increase between experimental groups to check if there are significant differences

# CF12
# Reshape data

CF12_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF12" & ( Day == "1" | Day == "15"))
CF12_incr_data$group <- c(rep("no_pOXA", 12), rep(c("pOXA", "pOXA+Ab"), 11))

CF12_incr_data <- CF12_incr_data[-34,] # delete outlier of failed curve

CF12_incr_data$Day <- as.factor(CF12_incr_data$Day)
CF12_incr_data$group <- as.factor(CF12_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = CF12_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF12_incr_data$AUC,CF12_incr_data$Day)

# Compute anova

CF12.aov <- aov(AUC ~ Day + group, data = CF12_incr_data); CF12.aov
summary(CF12.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(CF12_incr_data, aes(Day, AUC)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()

summary(glht(CF12.aov, linfct = mcp(group = "Tukey")))


# -------------------------------------------------------

## CF13

# Reshape data

CF13_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "CF13" & ( Day == "1" | Day == "15"))
CF13_incr_data$group <- c(rep("no_pOXA", 12), rep(c("pOXA", "pOXA+Ab"), 4), rep(c("pOXA+Ab","pOXA"),5), "pOXA+Ab", "pOXA+Ab")

CF13_incr_data <- CF13_incr_data[-34,] # delete outlier of failed curve

CF13_incr_data$Day <- as.factor(CF13_incr_data$Day)
CF13_incr_data$group <- as.factor(CF13_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = CF13_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(CF13_incr_data$AUC,CF13_incr_data$Day)

# Compute anova

CF13.wilk <- wilcox.test(AUC ~ Day, data = CF13_incr_data)
summary(CF13.wilk)

CF13.aov <- aov(AUC ~ Day, data = CF13_incr_data)
summary(CF13.aov)

# Plot data and make multiple comparisons by Tukey contrasts
ggplot(CF13_incr_data, aes(Day, AUC)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()

summary(glht(CF13.aov, linfct = mcp(group = "Tukey")))

# ----------------------------------------------------

## H53
# Reshape data

H53_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "H53" & ( Day == "1" | Day == "15"))
H53_incr_data$group <- c(rep("no_pOXA", 12),"pOXA","pOXA+Ab","pOXA+Ab","pOXA", rep(c("pOXA+Ab", "pOXA"), 4),"pOXA+Ab", rep(c("pOXA+Ab", "pOXA"), 4),"pOXA+Ab")

H53_incr_data$Day <- as.factor(H53_incr_data$Day)
H53_incr_data$group <- as.factor(H53_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = H53_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(H53_incr_data$AUC,H53_incr_data$Day)

# Compute anova

H53.aov <- aov(AUC ~ Day + group, data = H53_incr_data); H53.aov
summary(H53.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(H53_incr_data, aes(group, AUC, col = Day)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()

summary(glht(H53.aov, linfct = mcp(group = "Tukey")))

# ----------------------------------------------------------------

## K091
# Reshape data

K091_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "K091" & ( Day == "1" | Day == "15"))
K091_incr_data$group <- c(rep("no_pOXA", 12), rep(c("pOXA+Ab", "pOXA"), 2), rep(c("pOXA", "pOXA+Ab"), 3),"pOXA+Ab", "pOXA","pOXA+Ab", "pOXA", "pOXA", "pOXA+Ab", rep(c("pOXA", "pOXA+Ab"), 2))

K091_incr_data$Day <- as.factor(K091_incr_data$Day)
K091_incr_data$group <- as.factor(K091_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = K091_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K091_incr_data$AUC,K091_incr_data$Day)

# Compute anova

K091.aov <- aov(AUC ~ Day + group, data = K091_incr_data); K091.aov
summary(K091.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(K091_incr_data, aes(group, AUC, col = Day)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()
summary(glht(K091.aov, linfct = mcp(group = "Tukey")))

# -----------------------------------------------------------------------

# Reshape data

K147_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "K147" & ( Day == "1" | Day == "15"))
K147_incr_data$group <- c(rep("no_pOXA", 12), rep(c("pOXA", "pOXA+Ab"), 12))

K147_incr_data$Day <- as.factor(K147_incr_data$Day)
K147_incr_data$group <- as.factor(K147_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = K147_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K147_incr_data$AUC,K147_incr_data$Day)

# Compute anova

K147.aov <- aov(AUC ~ Day + group, data = K147_incr_data); K147.aov
summary(K147.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(K147_incr_data, aes(group, AUC, col = Day)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()
summary(glht(K147.aov, linfct = mcp(group = "Tukey")))

# ---------------------------------------------------------------

# Reshape data

K153_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "K153" & ( Day == "1" | Day == "15"))
K153_incr_data$group <- c(rep("no_pOXA", 12), rep(c("pOXA", "pOXA+Ab"), 11))

K153_incr_data <- K153_incr_data[-15,] # delete outlier of failed curve

K153_incr_data$Day <- as.factor(K153_incr_data$Day)
K153_incr_data$group <- as.factor(K153_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = K153_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K153_incr_data$AUC,K153_incr_data$Day)

# Compute anova

K153.aov <- aov(AUC ~ Day + group, data = K153_incr_data); K153.aov
summary(K153.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(K153_incr_data, aes(group, AUC, col = Day)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()
summary(glht(K153.aov, linfct = mcp(group = "Tukey")))

# ---------------------------------------------------------------

## K163
# Reshape data

K163_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "K163" & ( Day == "1" | Day == "15"))
K163_incr_data$group <- c(rep("no_pOXA", 12), rep(c("pOXA", "pOXA+Ab"), 11))

K163_incr_data <- K163_incr_data[-c(11,13,15,17),] # delete outlier of failed curve

K163_incr_data$Day <- as.factor(K163_incr_data$Day)
K163_incr_data$group <- as.factor(K163_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = K163_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K163_incr_data$AUC,K163_incr_data$Day)

# Compute anova

K163.aov <- aov(AUC ~ Day + group, data = K163_incr_data); K163.aov
summary(K163.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(K163_incr_data, aes(group, AUC, col = Day)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()
summary(glht(K163.aov, linfct = mcp(group = "Tukey")))

# --------------------------------------------------------------

## K209

# Reshape data

K209_incr_data <- whole_AUC_df_2 %>%
  filter(Sample == "K209" & ( Day == "1" | Day == "15"))
K209_incr_data$group <- c(rep("no_pOXA", 10), rep(c("pOXA", "pOXA+Ab"), 12))

K209_incr_data$Day <- as.factor(K209_incr_data$Day)
K209_incr_data$group <- as.factor(K209_incr_data$group)

# Check normality assumption

model  <- lm(AUC ~ Day + group, data = K209_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_incr_data$AUC,K209_incr_data$Day)

# Compute anova

K209.aov <- aov(AUC ~ Day + group, data = K209_incr_data); K209.aov
summary(K209.aov)

# Plot data and make multiple comparisons by Tukey contrasts

ggplot(K209_incr_data, aes(group, AUC, col = Day)) +
  geom_boxplot(aes(x = Day, y = AUC, col = group), width=1.5, alpha = 0.2, lwd = 1) +
  xlab("Day") + ylab("AUC") + theme(axis.title.y = element_text(size = 14)) + theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) + theme_bw()
summary(glht(K209.aov, linfct = mcp(group = "Tukey")))

# Now, for species

spc_incr_data <- rbind(CF12_incr_data, CF13_incr_data, H53_incr_data, K091_incr_data, K147_incr_data, K153_incr_data, K163_incr_data, K209_incr_data)
spc_incr_data <- spc_incr_data[-c(45,47,49,51,53,54,168,171,46,211, 50, 44),]

spc_incr_data$Species <- as.factor(spc_incr_data$Species)
spc_incr_data$Day <- as.factor(spc_incr_data$Day)


# Check normality assumption

model  <- lm(AUC ~ Species + Day, data = spc_incr_data)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

# Check homoscedasticity

leveneTest(K209_incr_data$AUC,K209_incr_data$Day)

# Compute anova
spc.wilc <- wilcox.test(AUC ~ Species, data = spc_incr_data)
summary(spc.wilc)
