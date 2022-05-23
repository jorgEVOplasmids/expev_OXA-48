setwd("C:/Users/jorge/OneDrive/Documentos/Universidad/TFM/clinical_samples/")

library(xlsx)
library(writexl)
library(sqldf)
library(dplyr)
library(data.table)

IS_rearrang <- read.xlsx(file="C:/Users/jorge/OneDrive/Documentos/Universidad/TFM/clinical_samples/IS_rearrangements.xlsx", sheetIndex=2,header = TRUE)

IS_rearrang$Strain <- NULL

shapiro.test(as.numeric(IS_rearrang[1,]))
shapiro.test(as.numeric(IS_rearrang[7,]))

t.test(IS_rearrang[7,], IS_rearrang[1,]) # pOXA-48 vs no pOXA-48
