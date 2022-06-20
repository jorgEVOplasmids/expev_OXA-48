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

path_to_txt= "txt/"
path_to_Growthrates="GrowthRates_results/"
path_to_output="Output/"
Rosett_Alfonso <- read.delim("Rosett/Rosett_Alfonso")

file.list <- list.files(path =path_to_txt, full.names = F)
df.list <- lapply(paste0(path_to_txt, file.list), 
                  function(x)read.delim(x, header=T, nrows=133, dec=","))

attr(df.list, "names") <- file.list
df <- bind_rows(df.list, .id = "id") %>% mutate(Time=rep(seq(0, (133*10)-10, 10),35))


df_test<-df %>% 
  select( -`Tâ.z.Optical.Density.600`) %>% gather(-Time,-id,  key = Well, value = OD ) %>% 
  separate(id, into=c("Plate", "Day"), remove = F) 

write.table(df_test, file=paste0(path_to_output, "curves_alf"))


curve_data <- df_test %>% 
  left_join(Rosett_Alfonso %>% mutate(Plate=as.character(Plate), 
                                      Day=as.character(Day),
                                      Well=as.character(Well))) %>% 
  mutate(Day=as.numeric(Day))

curve_data <- curve_data %>% 
  filter(
    Sample!="-") %>% 
  mutate(OD=as.numeric(OD))

#curvas RIBO
curve_data %>% 
  filter(Project=="RIBO") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  facet_grid(~Plasmid~Antibiotic~Day~Sample)+
  geom_hline(yintercept=0.6693333, linetype="dotted", color="darkgrey")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

test1 <- curve_data %>% filter(Sample=="Delta6",Day=="1",Time=="1300") %>% 
  group_by(Sample) %>% 
  mutate(mediaMaxOD=mean(OD))

curve_data %>% 
  filter(Sample=="K153") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K153")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))


curve_data %>% 
  filter(Sample=="K163") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K163")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K25") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K25")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C011") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C011")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C021") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C021")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C324") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C324")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C232") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C232")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C286") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C286")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="C309") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="C309")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K091") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K091")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K112") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K112")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K131") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K131")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K209") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K209")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="H53") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="H53")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="K147") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="K147")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="CF12") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="CF12")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Sample=="CF13") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day)+
  labs(title="CF13")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Project!="RIBO") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey")+
  facet_grid(~Plasmid~Antibiotic~Day~Sample)+
  labs(title="ALL")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

curve_data %>% 
  filter(Project=="RIBO") %>%
  ggplot(aes(y=OD, x=Time, color=Replicate, fill=Replicate, group=Replicate)) +
  geom_line()+
  facet_grid(~Plasmid~Antibiotic~Day~Sample)+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))


if (file.exists(paste0(path_to_Growthrates, "GrowthRates_results"))){
  Growthrate_results<-read.table((paste0(path_to_Growthrates, "GrowthRates_results")), header=T) %>% filter(r2>0.95)
  
} else{
  manysplits<- all_easylinear(OD~Time | Plasmid + Sample +  Antibiotic + Replicate + Day + Project + Species,
                              data=anti_join(curve_data, curve_data %>% filter(Sample=="CF13" & Day==13 & Antibiotic=="NO" & Plasmid=="NO")))
  write.table(results(manysplits), paste0(path_to_Growthrates, "GrowthRates_results"))
  
  Growthrate_results<-results(manysplits) %>% filter(r2>0.95)
}

data_analysed<- curve_data %>% group_by(Plasmid,Sample, Day,  Replicate, Antibiotic, Project, Species) %>% 
  group_modify(~ as.data.frame(flux::auc(.x$Time, .x$OD))) %>%
  mutate(AUC=`flux::auc(.x$Time, .x$OD)`) %>% 
  select(-`flux::auc(.x$Time, .x$OD)`)%>% 
  ungroup() 

data_analysed_2<- data_analysed %>% 
  left_join(curve_data %>% 
              group_by(Plasmid, Sample, Day, Antibiotic, Project, Replicate, Species) %>% 
              summarise(ODmax=max(OD, na.rm = T)))

data<-data_analysed_2 %>% 
  mutate(Replicate_cntrl=Replicate, 
         Replicate=as.numeric(Replicate)-1) %>% 
  left_join(Growthrate_results)

data<-data %>% 
  filter(Day==3) %>% #Este es el dÃ­a que usamos como referencia
  group_by(Plasmid, Sample, Day, Replicate, Antibiotic, Project, Species) %>% 
  summarise(Vmax_day3=mumax, 
            lag_day3=lag, 
            AUC_day3=AUC, 
            ODmax_day3=ODmax) %>% 
  left_join(data, by=c("Plasmid", "Sample", "Replicate", "Antibiotic", "Project", "Species")) %>% 
  select(-Day.x) %>% 
  mutate(Day=Day.y) %>% 
  mutate(Treatment=ifelse(Plasmid=="YES" & Antibiotic=="YES", "Plasmid+Ab",
                          ifelse(Plasmid=="YES" & Antibiotic=="NO", "Plasmid", 
                                 ifelse(Plasmid=="NO" & Antibiotic=="NO", "Control", "problems"))))

# Analysis of the delta in AUC comparing days 15 to 1 and 15 to 3

# Day 1

H53_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 1, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
H53_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 1, Plasmid == "YES", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
H53_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K091_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 4) %>%
  group_by(Sample)
K091_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 1, Plasmid == "YES", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K091_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 1, Plasmid == "YES", Antibiotic == "YES", Replicate != 3) %>%
  group_by(Sample)

K147_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)
K147_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K147_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K153_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 1, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K153_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K153_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K163_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
K163_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K163_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K209_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K209_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K209_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

C021_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
C021_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 1, Plasmid == "YES", Antibiotic == "NO", Replicate != 3) %>%
  group_by(Sample)
C021_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 1, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C324_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 4 | Replicate != 6) %>%
  group_by(Sample)
C324_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C324_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


C232_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 1, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C232_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C232_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 1, Plasmid == "YES", Antibiotic == "YES", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)


C286_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 2 | Replicate != 5) %>%
  group_by(Sample)
C286_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 1, Plasmid == "YES", Antibiotic == "NO", Replicate != 2 | Replicate != 4 | Replicate != 5) %>%
  group_by(Sample)
C286_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 1, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C309_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 1, Plasmid == "NO", Antibiotic == "NO",  Replicate != 1 | Replicate != 3 | Replicate != 5) %>%
  group_by(Sample)
C309_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C309_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 1, Plasmid == "YES", Antibiotic == "YES", Replicate != 4 | Replicate != 5 | Replicate != 6) %>%
  group_by(Sample)


K112_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 1, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K112_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 1, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K112_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF12_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 1, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
CF12_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 1, Plasmid == "YES", Antibiotic == "NO", Replicate != 5) %>%
  group_by(Sample)
CF12_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF13_no_no_day1 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 1, Plasmid == "NO", Antibiotic == "NO", Replicate != 4) %>%
  group_by(Sample)
CF13_no_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 1, Plasmid == "YES", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
CF13_yes_yes_day1 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 1, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

# Day 3

H53_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 3, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
H53_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 3, Plasmid == "YES", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
H53_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K091_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 4) %>%
  group_by(Sample)
K091_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 3, Plasmid == "YES", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K091_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 3, Plasmid == "YES", Antibiotic == "YES", Replicate != 3) %>%
  group_by(Sample)

K147_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)
K147_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K147_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K153_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 3, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K153_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K153_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K163_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
K163_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K163_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K209_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K209_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K209_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

C021_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
C021_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 3, Plasmid == "YES", Antibiotic == "NO", Replicate != 3) %>%
  group_by(Sample)
C021_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 3, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C324_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 4 | Replicate != 6) %>%
  group_by(Sample)
C324_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C324_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


C232_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 3, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C232_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C232_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 3, Plasmid == "YES", Antibiotic == "YES", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)


C286_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 2 | Replicate != 5) %>%
  group_by(Sample)
C286_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 3, Plasmid == "YES", Antibiotic == "NO", Replicate != 2 | Replicate != 4 | Replicate != 5) %>%
  group_by(Sample)
C286_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 3, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C309_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 3, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C309_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C309_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 3, Plasmid == "YES", Antibiotic == "YES", Replicate != 4 | Replicate != 5 | Replicate != 6) %>%
  group_by(Sample)


K112_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 3, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K112_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 3, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K112_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF12_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 3, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
CF12_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 3, Plasmid == "YES", Antibiotic == "NO", Replicate != 5) %>%
  group_by(Sample)
CF12_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF13_no_no_day3 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 3, Plasmid == "NO", Antibiotic == "NO", Replicate != 4) %>%
  group_by(Sample)
CF13_no_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 3, Plasmid == "YES", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
CF13_yes_yes_day3 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 3, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

# Day 5

H53_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 5, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
H53_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 5, Plasmid == "YES", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
H53_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K091_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 4) %>%
  group_by(Sample)
K091_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 5, Plasmid == "YES", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K091_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 5, Plasmid == "YES", Antibiotic == "YES", Replicate != 3) %>%
  group_by(Sample)

K147_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)
K147_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K147_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K153_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 5, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K153_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K153_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K163_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
K163_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K163_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K209_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K209_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K209_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

C021_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
C021_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 5, Plasmid == "YES", Antibiotic == "NO", Replicate != 3) %>%
  group_by(Sample)
C021_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 5, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C324_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 4 | Replicate != 6) %>%
  group_by(Sample)
C324_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C324_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


C232_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 5, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C232_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C232_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 5, Plasmid == "YES", Antibiotic == "YES", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)


C286_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 2 | Replicate != 5) %>%
  group_by(Sample)
C286_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 5, Plasmid == "YES", Antibiotic == "NO", Replicate != 2 | Replicate != 4 | Replicate != 5) %>%
  group_by(Sample)
C286_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 5, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C309_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 5, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C309_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C309_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 5, Plasmid == "YES", Antibiotic == "YES", Replicate != 4 | Replicate != 5 | Replicate != 6) %>%
  group_by(Sample)


K112_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 5, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K112_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 5, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K112_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF12_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 5, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
CF12_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 5, Plasmid == "YES", Antibiotic == "NO", Replicate != 5) %>%
  group_by(Sample)
CF12_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF13_no_no_day5 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 5, Plasmid == "NO", Antibiotic == "NO", Replicate != 4) %>%
  group_by(Sample)
CF13_no_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 5, Plasmid == "YES", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
CF13_yes_yes_day5 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 5, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

# Day 15

H53_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 15, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
H53_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 15, Plasmid == "YES", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
H53_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "H53", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K091_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 4) %>%
  group_by(Sample)
K091_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 15, Plasmid == "YES", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K091_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K091", Day == 15, Plasmid == "YES", Antibiotic == "YES", Replicate != 3) %>%
  group_by(Sample)

K147_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)
K147_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K147_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K147", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K153_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 15, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K153_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K153_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K153", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K163_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 2) %>%
  group_by(Sample)
K163_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K163_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K163", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

K209_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 1) %>%
  group_by(Sample)
K209_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K209_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K209", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

C021_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
C021_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 15, Plasmid == "YES", Antibiotic == "NO", Replicate != 3) %>%
  group_by(Sample)
C021_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C021", Day == 15, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C324_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 4 | Replicate != 6) %>%
  group_by(Sample)
C324_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C324_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C324", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


C232_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 15, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C232_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C232_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C232", Day == 15, Plasmid == "YES", Antibiotic == "YES", Replicate != 1 | Replicate != 2) %>%
  group_by(Sample)


C286_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 2 | Replicate != 5) %>%
  group_by(Sample)
C286_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 15, Plasmid == "YES", Antibiotic == "NO", Replicate != 2 | Replicate != 4 | Replicate != 5) %>%
  group_by(Sample)
C286_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C286", Day == 15, Plasmid == "YES", Antibiotic == "YES", Replicate != 5) %>%
  group_by(Sample)


C309_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 15, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
C309_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
C309_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "C309", Day == 15, Plasmid == "YES", Antibiotic == "YES", Replicate != 4 | Replicate != 5 | Replicate != 6) %>%
  group_by(Sample)


K112_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 15, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
K112_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 15, Plasmid == "YES", Antibiotic == "NO") %>%
  group_by(Sample)
K112_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "K112", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF12_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 15, Plasmid == "NO", Antibiotic == "NO") %>%
  group_by(Sample)
CF12_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 15, Plasmid == "YES", Antibiotic == "NO", Replicate != 5) %>%
  group_by(Sample)
CF12_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "CF12", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)


CF13_no_no_day15 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 15, Plasmid == "NO", Antibiotic == "NO", Replicate != 4) %>%
  group_by(Sample)
CF13_no_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 15, Plasmid == "YES", Antibiotic == "NO", Replicate != 6) %>%
  group_by(Sample)
CF13_yes_yes_day15 <- data_analysed_2 %>%
  filter(Sample == "CF13", Day == 15, Plasmid == "YES", Antibiotic == "YES") %>%
  group_by(Sample)

# Calculate the mean values of samples without pOXA-48 per strain to use them as reference to calculate the costs, and also the median values

H53_no_no_day1$meanAUC <- mean(H53_no_no_day1$AUC)
K091_no_no_day1$meanAUC <- mean(K091_no_no_day1$AUC)
K147_no_no_day1$meanAUC <- mean(K147_no_no_day1$AUC)
K153_no_no_day1$meanAUC <- mean(K153_no_no_day1$AUC)
K163_no_no_day1$meanAUC <- mean(K163_no_no_day1$AUC)
K209_no_no_day1$meanAUC <- mean(K209_no_no_day1$AUC)

H53_no_no_day1$medianAUC <- median(H53_no_no_day1$AUC)
K091_no_no_day1$medianAUC <- median(K091_no_no_day1$AUC)
K147_no_no_day1$medianAUC <- median(K147_no_no_day1$AUC)
K153_no_no_day1$medianAUC <- median(K153_no_no_day1$AUC)
K163_no_no_day1$medianAUC <- median(K163_no_no_day1$AUC)
K209_no_no_day1$medianAUC <- median(K209_no_no_day1$AUC)

# And also the mean and median values of AUC in samples in days 3, 5 and 15

H53_no_no_day3$meanAUC <- mean(H53_no_no_day3$AUC)
K091_no_no_day3$meanAUC <- mean(K091_no_no_day3$AUC)
K147_no_no_day3$meanAUC <- mean(K147_no_no_day3$AUC)
K153_no_no_day3$meanAUC <- mean(K153_no_no_day3$AUC)
K163_no_no_day3$meanAUC <- mean(K163_no_no_day3$AUC)
K209_no_no_day3$meanAUC <- mean(K209_no_no_day3$AUC)

H53_no_no_day3$medianAUC <- median(H53_no_no_day3$AUC)
K091_no_no_day3$medianAUC <- median(K091_no_no_day3$AUC)
K147_no_no_day3$medianAUC <- median(K147_no_no_day3$AUC)
K153_no_no_day3$medianAUC <- median(K153_no_no_day3$AUC)
K163_no_no_day3$medianAUC <- median(K163_no_no_day3$AUC)
K209_no_no_day3$medianAUC <- median(K209_no_no_day3$AUC)

H53_no_no_day5$meanAUC <- mean(H53_no_no_day5$AUC)
K091_no_no_day5$meanAUC <- mean(K091_no_no_day5$AUC)
K147_no_no_day5$meanAUC <- mean(K147_no_no_day5$AUC)
K153_no_no_day5$meanAUC <- mean(K153_no_no_day5$AUC)
K163_no_no_day5$meanAUC <- mean(K163_no_no_day5$AUC)
K209_no_no_day5$meanAUC <- mean(K209_no_no_day5$AUC)

H53_no_no_day5$medianAUC <- median(H53_no_no_day5$AUC)
K091_no_no_day5$medianAUC <- median(K091_no_no_day5$AUC)
K147_no_no_day5$medianAUC <- median(K147_no_no_day5$AUC)
K153_no_no_day5$medianAUC <- median(K153_no_no_day5$AUC)
K163_no_no_day5$medianAUC <- median(K163_no_no_day5$AUC)
K209_no_no_day5$medianAUC <- median(K209_no_no_day5$AUC)

H53_no_no_day15$meanAUC <- mean(H53_no_no_day15$AUC)
K091_no_no_day15$meanAUC <- mean(K091_no_no_day15$AUC)
K147_no_no_day15$meanAUC <- mean(K147_no_no_day15$AUC)
K153_no_no_day15$meanAUC <- mean(K153_no_no_day15$AUC)
K163_no_no_day15$meanAUC <- mean(K163_no_no_day15$AUC)
K209_no_no_day15$meanAUC <- mean(K209_no_no_day15$AUC)

H53_no_no_day15$medianAUC <- median(H53_no_no_day15$AUC)
K091_no_no_day15$medianAUC <- median(K091_no_no_day15$AUC)
K147_no_no_day15$medianAUC <- median(K147_no_no_day15$AUC)
K153_no_no_day15$medianAUC <- median(K153_no_no_day15$AUC)
K163_no_no_day15$medianAUC <- median(K163_no_no_day15$AUC)
K209_no_no_day15$medianAUC <- median(K209_no_no_day15$AUC)

# That has to be done for the other 2 conditions

H53_no_yes_day1$meanAUC <- mean(H53_no_yes_day1$AUC)
K091_no_yes_day1$meanAUC <- mean(K091_no_yes_day1$AUC)
K147_no_yes_day1$meanAUC <- mean(K147_no_yes_day1$AUC)
K153_no_yes_day1$meanAUC <- mean(K153_no_yes_day1$AUC)
K163_no_yes_day1$meanAUC <- mean(K163_no_yes_day1$AUC)
K209_no_yes_day1$meanAUC <- mean(K209_no_yes_day1$AUC)

H53_no_yes_day1$medianAUC <- median(H53_no_yes_day1$AUC)
K091_no_yes_day1$medianAUC <- median(K091_no_yes_day1$AUC)
K147_no_yes_day1$medianAUC <- median(K147_no_yes_day1$AUC)
K153_no_yes_day1$medianAUC <- median(K153_no_yes_day1$AUC)
K163_no_yes_day1$medianAUC <- median(K163_no_yes_day1$AUC)
K209_no_yes_day1$medianAUC <- median(K209_no_yes_day1$AUC)

H53_no_yes_day3$meanAUC <- mean(H53_no_yes_day3$AUC)
K091_no_yes_day3$meanAUC <- mean(K091_no_yes_day3$AUC)
K147_no_yes_day3$meanAUC <- mean(K147_no_yes_day3$AUC)
K153_no_yes_day3$meanAUC <- mean(K153_no_yes_day3$AUC)
K163_no_yes_day3$meanAUC <- mean(K163_no_yes_day3$AUC)
K209_no_yes_day3$meanAUC <- mean(K209_no_yes_day3$AUC)

H53_no_yes_day3$medianAUC <- median(H53_no_yes_day3$AUC)
K091_no_yes_day3$medianAUC <- median(K091_no_yes_day3$AUC)
K147_no_yes_day3$medianAUC <- median(K147_no_yes_day3$AUC)
K153_no_yes_day3$medianAUC <- median(K153_no_yes_day3$AUC)
K163_no_yes_day3$medianAUC <- median(K163_no_yes_day3$AUC)
K209_no_yes_day3$medianAUC <- median(K209_no_yes_day3$AUC)

H53_no_yes_day5$meanAUC <- mean(H53_no_yes_day5$AUC)
K091_no_yes_day5$meanAUC <- mean(K091_no_yes_day5$AUC)
K147_no_yes_day5$meanAUC <- mean(K147_no_yes_day5$AUC)
K153_no_yes_day5$meanAUC <- mean(K153_no_yes_day5$AUC)
K163_no_yes_day5$meanAUC <- mean(K163_no_yes_day5$AUC)
K209_no_yes_day5$meanAUC <- mean(K209_no_yes_day5$AUC)

H53_no_yes_day5$medianAUC <- median(H53_no_yes_day5$AUC)
K091_no_yes_day5$medianAUC <- median(K091_no_yes_day5$AUC)
K147_no_yes_day5$medianAUC <- median(K147_no_yes_day5$AUC)
K153_no_yes_day5$medianAUC <- median(K153_no_yes_day5$AUC)
K163_no_yes_day5$medianAUC <- median(K163_no_yes_day5$AUC)
K209_no_yes_day5$medianAUC <- median(K209_no_yes_day5$AUC)

H53_no_yes_day15$meanAUC <- mean(H53_no_yes_day15$AUC)
K091_no_yes_day15$meanAUC <- mean(K091_no_yes_day15$AUC)
K147_no_yes_day15$meanAUC <- mean(K147_no_yes_day15$AUC)
K153_no_yes_day15$meanAUC <- mean(K153_no_yes_day15$AUC)
K163_no_yes_day15$meanAUC <- mean(K163_no_yes_day15$AUC)
K209_no_yes_day15$meanAUC <- mean(K209_no_yes_day15$AUC)

H53_no_yes_day15$medianAUC <- median(H53_no_yes_day15$AUC)
K091_no_yes_day15$medianAUC <- median(K091_no_yes_day15$AUC)
K147_no_yes_day15$medianAUC <- median(K147_no_yes_day15$AUC)
K153_no_yes_day15$medianAUC <- median(K153_no_yes_day15$AUC)
K163_no_yes_day15$medianAUC <- median(K163_no_yes_day15$AUC)
K209_no_yes_day15$medianAUC <- median(K209_no_yes_day15$AUC)

H53_yes_yes_day1$meanAUC <- mean(H53_yes_yes_day1$AUC)
K091_yes_yes_day1$meanAUC <- mean(K091_yes_yes_day1$AUC)
K147_yes_yes_day1$meanAUC <- mean(K147_yes_yes_day1$AUC)
K153_yes_yes_day1$meanAUC <- mean(K153_yes_yes_day1$AUC)
K163_yes_yes_day1$meanAUC <- mean(K163_yes_yes_day1$AUC)
K209_yes_yes_day1$meanAUC <- mean(K209_yes_yes_day1$AUC)

H53_yes_yes_day1$medianAUC <- median(H53_yes_yes_day1$AUC)
K091_yes_yes_day1$medianAUC <- median(K091_yes_yes_day1$AUC)
K147_yes_yes_day1$medianAUC <- median(K147_yes_yes_day1$AUC)
K153_yes_yes_day1$medianAUC <- median(K153_yes_yes_day1$AUC)
K163_yes_yes_day1$medianAUC <- median(K163_yes_yes_day1$AUC)
K209_yes_yes_day1$medianAUC <- median(K209_yes_yes_day1$AUC)

H53_yes_yes_day3$meanAUC <- mean(H53_yes_yes_day3$AUC)
K091_yes_yes_day3$meanAUC <- mean(K091_yes_yes_day3$AUC)
K147_yes_yes_day3$meanAUC <- mean(K147_yes_yes_day3$AUC)
K153_yes_yes_day3$meanAUC <- mean(K153_yes_yes_day3$AUC)
K163_yes_yes_day3$meanAUC <- mean(K163_yes_yes_day3$AUC)
K209_yes_yes_day3$meanAUC <- mean(K209_yes_yes_day3$AUC)

H53_yes_yes_day3$medianAUC <- median(H53_yes_yes_day3$AUC)
K091_yes_yes_day3$medianAUC <- median(K091_yes_yes_day3$AUC)
K147_yes_yes_day3$medianAUC <- median(K147_yes_yes_day3$AUC)
K153_yes_yes_day3$medianAUC <- median(K153_yes_yes_day3$AUC)
K163_yes_yes_day3$medianAUC <- median(K163_yes_yes_day3$AUC)
K209_yes_yes_day3$medianAUC <- median(K209_yes_yes_day3$AUC)

H53_yes_yes_day5$meanAUC <- mean(H53_yes_yes_day5$AUC)
K091_yes_yes_day5$meanAUC <- mean(K091_yes_yes_day5$AUC)
K147_yes_yes_day5$meanAUC <- mean(K147_yes_yes_day5$AUC)
K153_yes_yes_day5$meanAUC <- mean(K153_yes_yes_day5$AUC)
K163_yes_yes_day5$meanAUC <- mean(K163_yes_yes_day5$AUC)
K209_yes_yes_day5$meanAUC <- mean(K209_yes_yes_day5$AUC)

H53_yes_yes_day5$medianAUC <- median(H53_yes_yes_day5$AUC)
K091_yes_yes_day5$medianAUC <- median(K091_yes_yes_day5$AUC)
K147_yes_yes_day5$medianAUC <- median(K147_yes_yes_day5$AUC)
K153_yes_yes_day5$medianAUC <- median(K153_yes_yes_day5$AUC)
K163_yes_yes_day5$medianAUC <- median(K163_yes_yes_day5$AUC)
K209_yes_yes_day5$medianAUC <- median(K209_yes_yes_day5$AUC)

H53_yes_yes_day15$meanAUC <- mean(H53_yes_yes_day15$AUC)
K091_yes_yes_day15$meanAUC <- mean(K091_yes_yes_day15$AUC)
K147_yes_yes_day15$meanAUC <- mean(K147_yes_yes_day15$AUC)
K153_yes_yes_day15$meanAUC <- mean(K153_yes_yes_day15$AUC)
K163_yes_yes_day15$meanAUC <- mean(K163_yes_yes_day15$AUC)
K209_yes_yes_day15$meanAUC <- mean(K209_yes_yes_day15$AUC)

H53_yes_yes_day15$medianAUC <- median(H53_yes_yes_day15$AUC)
K091_yes_yes_day15$medianAUC <- median(K091_yes_yes_day15$AUC)
K147_yes_yes_day15$medianAUC <- median(K147_yes_yes_day15$AUC)
K153_yes_yes_day15$medianAUC <- median(K153_yes_yes_day15$AUC)
K163_yes_yes_day15$medianAUC <- median(K163_yes_yes_day15$AUC)
K209_yes_yes_day15$medianAUC <- median(K209_yes_yes_day15$AUC)

# For the rest of the samples

C021_no_no_day1$meanAUC <- mean(C021_no_no_day1$AUC)
C324_no_no_day1$meanAUC <- mean(C324_no_no_day1$AUC)
C232_no_no_day1$meanAUC <- mean(C232_no_no_day1$AUC)
C286_no_no_day1$meanAUC <- mean(C286_no_no_day1$AUC)
C309_no_no_day1$meanAUC <- 816.2 # due to fail in calculous, manually introduced
K112_no_no_day1$meanAUC <- mean(K112_no_no_day1$AUC)
CF12_no_no_day1$meanAUC <- mean(CF12_no_no_day1$AUC)
CF13_no_no_day1$meanAUC <- mean(CF13_no_no_day1$AUC)

C021_no_no_day1$medianAUC <- median(C021_no_no_day1$AUC)
C324_no_no_day1$medianAUC <- median(C324_no_no_day1$AUC)
C232_no_no_day1$medianAUC <- median(C232_no_no_day1$AUC)
C286_no_no_day1$medianAUC <- median(C286_no_no_day1$AUC)
C309_no_no_day1$medianAUC <- median(C309_no_no_day1$AUC)
K112_no_no_day1$medianAUC <- median(K112_no_no_day1$AUC)
CF12_no_no_day1$medianAUC <- median(CF12_no_no_day1$AUC)
CF13_no_no_day1$medianAUC <- median(CF13_no_no_day1$AUC)

# And also the mean and median values of AUC in samples in days 3, 5 and 15

C021_no_no_day3$meanAUC <- mean(C021_no_no_day3$AUC)
C324_no_no_day3$meanAUC <- mean(C324_no_no_day3$AUC)
C232_no_no_day3$meanAUC <- mean(C232_no_no_day3$AUC)
C286_no_no_day3$meanAUC <- mean(C286_no_no_day3$AUC)
C309_no_no_day3$meanAUC <- mean(C309_no_no_day3$AUC)
K112_no_no_day3$meanAUC <- mean(K112_no_no_day3$AUC)
CF12_no_no_day3$meanAUC <- mean(CF12_no_no_day3$AUC)
CF13_no_no_day3$meanAUC <- mean(CF13_no_no_day3$AUC)

C021_no_no_day3$medianAUC <- median(C021_no_no_day3$AUC)
C324_no_no_day3$medianAUC <- median(C324_no_no_day3$AUC)
C232_no_no_day3$medianAUC <- median(C232_no_no_day3$AUC)
C286_no_no_day3$medianAUC <- median(C286_no_no_day3$AUC)
C309_no_no_day3$medianAUC <- median(C309_no_no_day3$AUC)
K112_no_no_day3$medianAUC <- median(K112_no_no_day3$AUC)
CF12_no_no_day3$medianAUC <- median(CF12_no_no_day3$AUC)
CF13_no_no_day3$medianAUC <- median(CF13_no_no_day3$AUC)

C021_no_no_day5$meanAUC <- mean(C021_no_no_day5$AUC)
C324_no_no_day5$meanAUC <- mean(C324_no_no_day5$AUC)
C232_no_no_day5$meanAUC <- mean(C232_no_no_day5$AUC)
C286_no_no_day5$meanAUC <- mean(C286_no_no_day5$AUC)
C309_no_no_day5$meanAUC <- mean(C309_no_no_day5$AUC)
K112_no_no_day5$meanAUC <- mean(K112_no_no_day5$AUC)
CF12_no_no_day5$meanAUC <- mean(CF12_no_no_day5$AUC)
CF13_no_no_day5$meanAUC <- mean(CF13_no_no_day5$AUC)

C021_no_no_day5$medianAUC <- median(C021_no_no_day5$AUC)
C324_no_no_day5$medianAUC <- median(C324_no_no_day5$AUC)
C232_no_no_day5$medianAUC <- median(C232_no_no_day5$AUC)
C286_no_no_day5$medianAUC <- median(C286_no_no_day5$AUC)
C309_no_no_day5$medianAUC <- median(C309_no_no_day5$AUC)
K112_no_no_day5$medianAUC <- median(K112_no_no_day5$AUC)
CF12_no_no_day5$medianAUC <- median(CF12_no_no_day5$AUC)
CF13_no_no_day5$medianAUC <- median(CF13_no_no_day5$AUC)

C021_no_no_day15$meanAUC <- mean(C021_no_no_day15$AUC)
C324_no_no_day15$meanAUC <- mean(C324_no_no_day15$AUC)
C232_no_no_day15$meanAUC <- mean(C232_no_no_day15$AUC)
C286_no_no_day15$meanAUC <- mean(C286_no_no_day15$AUC)
C309_no_no_day15$meanAUC <- mean(C309_no_no_day15$AUC)
K112_no_no_day15$meanAUC <- mean(K112_no_no_day15$AUC)
CF12_no_no_day15$meanAUC <- mean(CF12_no_no_day15$AUC)
CF13_no_no_day15$meanAUC <- mean(CF13_no_no_day15$AUC)

C021_no_no_day15$medianAUC <- median(C021_no_no_day15$AUC)
C324_no_no_day15$medianAUC <- median(C324_no_no_day15$AUC)
C232_no_no_day15$medianAUC <- median(C232_no_no_day15$AUC)
C286_no_no_day15$medianAUC <- median(C286_no_no_day15$AUC)
C309_no_no_day15$medianAUC <- median(C309_no_no_day15$AUC)
K112_no_no_day15$medianAUC <- median(K112_no_no_day15$AUC)
CF12_no_no_day15$medianAUC <- median(CF12_no_no_day15$AUC)
CF13_no_no_day15$medianAUC <- median(CF13_no_no_day15$AUC)

# That has to be done for the other 2 conditions

C021_no_yes_day1$meanAUC <- mean(C021_no_yes_day1$AUC)
C324_no_yes_day1$meanAUC <- mean(C324_no_yes_day1$AUC)
C232_no_yes_day1$meanAUC <- mean(C232_no_yes_day1$AUC)
C286_no_yes_day1$meanAUC <- mean(C286_no_yes_day1$AUC)
C309_no_yes_day1$meanAUC <- mean(C309_no_yes_day1$AUC)
K112_no_yes_day1$meanAUC <- mean(K112_no_yes_day1$AUC)
CF12_no_yes_day1$meanAUC <- mean(CF12_no_yes_day1$AUC)
CF13_no_yes_day1$meanAUC <- mean(CF13_no_yes_day1$AUC)

C021_no_yes_day1$medianAUC <- median(C021_no_yes_day1$AUC)
C324_no_yes_day1$medianAUC <- median(C324_no_yes_day1$AUC)
C232_no_yes_day1$medianAUC <- median(C232_no_yes_day1$AUC)
C286_no_yes_day1$medianAUC <- median(C286_no_yes_day1$AUC)
C309_no_yes_day1$medianAUC <- median(C309_no_yes_day1$AUC)
K112_no_yes_day1$medianAUC <- median(K112_no_yes_day1$AUC)
CF12_no_yes_day1$medianAUC <- median(CF12_no_yes_day1$AUC)
CF13_no_yes_day1$medianAUC <- median(CF13_no_yes_day1$AUC)

# And also the mean and median values of AUC in samples in days 3, 5 and 15

C021_no_yes_day3$meanAUC <- mean(C021_no_yes_day3$AUC)
C324_no_yes_day3$meanAUC <- mean(C324_no_yes_day3$AUC)
C232_no_yes_day3$meanAUC <- mean(C232_no_yes_day3$AUC)
C286_no_yes_day3$meanAUC <- mean(C286_no_yes_day3$AUC)
C309_no_yes_day3$meanAUC <- mean(C309_no_yes_day3$AUC)
K112_no_yes_day3$meanAUC <- mean(K112_no_yes_day3$AUC)
CF12_no_yes_day3$meanAUC <- mean(CF12_no_yes_day3$AUC)
CF13_no_yes_day3$meanAUC <- mean(CF13_no_yes_day3$AUC)

C021_no_yes_day3$medianAUC <- median(C021_no_yes_day3$AUC)
C324_no_yes_day3$medianAUC <- median(C324_no_yes_day3$AUC)
C232_no_yes_day3$medianAUC <- median(C232_no_yes_day3$AUC)
C286_no_yes_day3$medianAUC <- median(C286_no_yes_day3$AUC)
C309_no_yes_day3$medianAUC <- median(C309_no_yes_day3$AUC)
K112_no_yes_day3$medianAUC <- median(K112_no_yes_day3$AUC)
CF12_no_yes_day3$medianAUC <- median(CF12_no_yes_day3$AUC)
CF13_no_yes_day3$medianAUC <- median(CF13_no_yes_day3$AUC)

C021_no_yes_day5$meanAUC <- mean(C021_no_yes_day5$AUC)
C324_no_yes_day5$meanAUC <- mean(C324_no_yes_day5$AUC)
C232_no_yes_day5$meanAUC <- mean(C232_no_yes_day5$AUC)
C286_no_yes_day5$meanAUC <- mean(C286_no_yes_day5$AUC)
C309_no_yes_day5$meanAUC <- mean(C309_no_yes_day5$AUC)
K112_no_yes_day5$meanAUC <- mean(K112_no_yes_day5$AUC)
CF12_no_yes_day5$meanAUC <- mean(CF12_no_yes_day5$AUC)
CF13_no_yes_day5$meanAUC <- mean(CF13_no_yes_day5$AUC)

C021_no_yes_day5$medianAUC <- median(C021_no_yes_day5$AUC)
C324_no_yes_day5$medianAUC <- median(C324_no_yes_day5$AUC)
C232_no_yes_day5$medianAUC <- median(C232_no_yes_day5$AUC)
C286_no_yes_day5$medianAUC <- median(C286_no_yes_day5$AUC)
C309_no_yes_day5$medianAUC <- median(C309_no_yes_day5$AUC)
K112_no_yes_day5$medianAUC <- median(K112_no_yes_day5$AUC)
CF12_no_yes_day5$medianAUC <- median(CF12_no_yes_day5$AUC)
CF13_no_yes_day5$medianAUC <- median(CF13_no_yes_day5$AUC)

C021_no_yes_day15$meanAUC <- mean(C021_no_yes_day15$AUC)
C324_no_yes_day15$meanAUC <- mean(C324_no_yes_day15$AUC)
C232_no_yes_day15$meanAUC <- mean(C232_no_yes_day15$AUC)
C286_no_yes_day15$meanAUC <- mean(C286_no_yes_day15$AUC)
C309_no_yes_day15$meanAUC <- mean(C309_no_yes_day15$AUC)
K112_no_yes_day15$meanAUC <- mean(K112_no_yes_day15$AUC)
CF12_no_yes_day15$meanAUC <- mean(CF12_no_yes_day15$AUC)
CF13_no_yes_day15$meanAUC <- mean(CF13_no_yes_day15$AUC)

C021_no_yes_day15$medianAUC <- median(C021_no_yes_day15$AUC)
C324_no_yes_day15$medianAUC <- median(C324_no_yes_day15$AUC)
C232_no_yes_day15$medianAUC <- median(C232_no_yes_day15$AUC)
C286_no_yes_day15$medianAUC <- median(C286_no_yes_day15$AUC)
C309_no_yes_day15$medianAUC <- median(C309_no_yes_day15$AUC)
K112_no_yes_day15$medianAUC <- median(K112_no_yes_day15$AUC)
CF12_no_yes_day15$medianAUC <- median(CF12_no_yes_day15$AUC)
CF13_no_yes_day15$medianAUC <- median(CF13_no_yes_day15$AUC)

C021_yes_yes_day1$meanAUC <- mean(C021_yes_yes_day1$AUC)
C324_yes_yes_day1$meanAUC <- mean(C324_yes_yes_day1$AUC)
C232_yes_yes_day1$meanAUC <- mean(C232_yes_yes_day1$AUC)
C286_yes_yes_day1$meanAUC <- mean(C286_yes_yes_day1$AUC)
C309_yes_yes_day1$meanAUC <- mean(C309_yes_yes_day1$AUC)
K112_yes_yes_day1$meanAUC <- mean(K112_yes_yes_day1$AUC)
CF12_yes_yes_day1$meanAUC <- mean(CF12_yes_yes_day1$AUC)
CF13_yes_yes_day1$meanAUC <- mean(CF13_yes_yes_day1$AUC)

C021_yes_yes_day1$medianAUC <- median(C021_yes_yes_day1$AUC)
C324_yes_yes_day1$medianAUC <- median(C324_yes_yes_day1$AUC)
C232_yes_yes_day1$medianAUC <- median(C232_yes_yes_day1$AUC)
C286_yes_yes_day1$medianAUC <- median(C286_yes_yes_day1$AUC)
C309_yes_yes_day1$medianAUC <- median(C309_yes_yes_day1$AUC)
K112_yes_yes_day1$medianAUC <- median(K112_yes_yes_day1$AUC)
CF12_yes_yes_day1$medianAUC <- median(CF12_yes_yes_day1$AUC)
CF13_yes_yes_day1$medianAUC <- median(CF13_yes_yes_day1$AUC)

# And also the mean and median values of AUC in samples in days 3, 5 and 15

C021_yes_yes_day3$meanAUC <- mean(C021_yes_yes_day3$AUC)
C324_yes_yes_day3$meanAUC <- mean(C324_yes_yes_day3$AUC)
C232_yes_yes_day3$meanAUC <- mean(C232_yes_yes_day3$AUC)
C286_yes_yes_day3$meanAUC <- mean(C286_yes_yes_day3$AUC)
C309_yes_yes_day3$meanAUC <- mean(C309_yes_yes_day3$AUC)
K112_yes_yes_day3$meanAUC <- mean(K112_yes_yes_day3$AUC)
CF12_yes_yes_day3$meanAUC <- mean(CF12_yes_yes_day3$AUC)
CF13_yes_yes_day3$meanAUC <- mean(CF13_yes_yes_day3$AUC)

C021_yes_yes_day3$medianAUC <- median(C021_yes_yes_day3$AUC)
C324_yes_yes_day3$medianAUC <- median(C324_yes_yes_day3$AUC)
C232_yes_yes_day3$medianAUC <- median(C232_yes_yes_day3$AUC)
C286_yes_yes_day3$medianAUC <- median(C286_yes_yes_day3$AUC)
C309_yes_yes_day3$medianAUC <- median(C309_yes_yes_day3$AUC)
K112_yes_yes_day3$medianAUC <- median(K112_yes_yes_day3$AUC)
CF12_yes_yes_day3$medianAUC <- median(CF12_yes_yes_day3$AUC)
CF13_yes_yes_day3$medianAUC <- median(CF13_yes_yes_day3$AUC)

C021_yes_yes_day5$meanAUC <- mean(C021_yes_yes_day5$AUC)
C324_yes_yes_day5$meanAUC <- mean(C324_yes_yes_day5$AUC)
C232_yes_yes_day5$meanAUC <- mean(C232_yes_yes_day5$AUC)
C286_yes_yes_day5$meanAUC <- mean(C286_yes_yes_day5$AUC)
C309_yes_yes_day5$meanAUC <- mean(C309_yes_yes_day5$AUC)
K112_yes_yes_day5$meanAUC <- mean(K112_yes_yes_day5$AUC)
CF12_yes_yes_day5$meanAUC <- mean(CF12_yes_yes_day5$AUC)
CF13_yes_yes_day5$meanAUC <- mean(CF13_yes_yes_day5$AUC)

C021_yes_yes_day5$medianAUC <- median(C021_yes_yes_day5$AUC)
C324_yes_yes_day5$medianAUC <- median(C324_yes_yes_day5$AUC)
C232_yes_yes_day5$medianAUC <- median(C232_yes_yes_day5$AUC)
C286_yes_yes_day5$medianAUC <- median(C286_yes_yes_day5$AUC)
C309_yes_yes_day5$medianAUC <- median(C309_yes_yes_day5$AUC)
K112_yes_yes_day5$medianAUC <- median(K112_yes_yes_day5$AUC)
CF12_yes_yes_day5$medianAUC <- median(CF12_yes_yes_day5$AUC)
CF13_yes_yes_day5$medianAUC <- median(CF13_yes_yes_day5$AUC)

C021_yes_yes_day15$meanAUC <- mean(C021_yes_yes_day15$AUC)
C324_yes_yes_day15$meanAUC <- mean(C324_yes_yes_day15$AUC)
C232_yes_yes_day15$meanAUC <- mean(C232_yes_yes_day15$AUC)
C286_yes_yes_day15$meanAUC <- mean(C286_yes_yes_day15$AUC)
C309_yes_yes_day15$meanAUC <- mean(C309_yes_yes_day15$AUC)
K112_yes_yes_day15$meanAUC <- mean(K112_yes_yes_day15$AUC)
CF12_yes_yes_day15$meanAUC <- mean(CF12_yes_yes_day15$AUC)
CF13_yes_yes_day15$meanAUC <- mean(CF13_yes_yes_day15$AUC)

C021_yes_yes_day15$medianAUC <- median(C021_yes_yes_day15$AUC)
C324_yes_yes_day15$medianAUC <- median(C324_yes_yes_day15$AUC)
C232_yes_yes_day15$medianAUC <- median(C232_yes_yes_day15$AUC)
C286_yes_yes_day15$medianAUC <- median(C286_yes_yes_day15$AUC)
C309_yes_yes_day15$medianAUC <- median(C309_yes_yes_day15$AUC)
K112_yes_yes_day15$medianAUC <- median(K112_yes_yes_day15$AUC)
CF12_yes_yes_day15$medianAUC <- median(CF12_yes_yes_day15$AUC)
CF13_yes_yes_day15$medianAUC <- median(CF13_yes_yes_day15$AUC)

# Lastly, include a column with the difference between the AUC of each replicate with/without Ab and with pOXA-48 and the mean of the strain without pOXA-48

H53_no_no_day1$relative_AUC <- (H53_no_no_day1$AUC)/H53_no_no_day1$meanAUC
K091_no_no_day1$relative_AUC <- (K091_no_no_day1$AUC)/K091_no_no_day1$meanAUC
K147_no_no_day1$relative_AUC <- (K147_no_no_day1$AUC)/K147_no_no_day1$meanAUC
K153_no_no_day1$relative_AUC <- (K153_no_no_day1$AUC)/K153_no_no_day1$meanAUC
K163_no_no_day1$relative_AUC <- (K163_no_no_day1$AUC)/K163_no_no_day1$meanAUC
K209_no_no_day1$relative_AUC <- (K209_no_no_day1$AUC)/K209_no_no_day1$meanAUC

H53_no_no_day3$relative_AUC <- (H53_no_no_day3$AUC)/H53_no_no_day1$meanAUC
K091_no_no_day3$relative_AUC <- (K091_no_no_day3$AUC)/K091_no_no_day1$meanAUC
K147_no_no_day3$relative_AUC <- (K147_no_no_day3$AUC)/K147_no_no_day1$meanAUC
K153_no_no_day3$relative_AUC <- (K153_no_no_day3$AUC)/K153_no_no_day1$meanAUC
K163_no_no_day3$relative_AUC <- (K163_no_no_day3$AUC)/K163_no_no_day1$meanAUC
K209_no_no_day3$relative_AUC <- (K209_no_no_day3$AUC)/K209_no_no_day1$meanAUC

H53_no_no_day5$relative_AUC <- (H53_no_no_day5$AUC)/H53_no_no_day1$meanAUC
K091_no_no_day5$relative_AUC <- (K091_no_no_day5$AUC)/K091_no_no_day1$meanAUC
K147_no_no_day5$relative_AUC <- (K147_no_no_day5$AUC)/K147_no_no_day1$meanAUC
K153_no_no_day5$relative_AUC <- (K153_no_no_day5$AUC)/K153_no_no_day1$meanAUC
K163_no_no_day5$relative_AUC <- (K163_no_no_day5$AUC)/K163_no_no_day1$meanAUC
K209_no_no_day5$relative_AUC <- (K209_no_no_day5$AUC)/K209_no_no_day1$meanAUC

H53_no_no_day15$relative_AUC <- (H53_no_no_day15$AUC)/H53_no_no_day1$meanAUC
K091_no_no_day15$relative_AUC <- (K091_no_no_day15$AUC)/K091_no_no_day1$meanAUC
K147_no_no_day15$relative_AUC <- (K147_no_no_day15$AUC)/K147_no_no_day1$meanAUC
K153_no_no_day15$relative_AUC <- (K153_no_no_day15$AUC)/K153_no_no_day1$meanAUC
K163_no_no_day15$relative_AUC <- (K163_no_no_day15$AUC)/K163_no_no_day1$meanAUC
K209_no_no_day15$relative_AUC <- (K209_no_no_day15$AUC)/K209_no_no_day1$meanAUC

# For samples with pOXA-48 in presence/absence of antibiotic

H53_no_yes_day1$relative_AUC <- (H53_no_yes_day1$AUC)/H53_no_no_day1$meanAUC[1]
K091_no_yes_day1$relative_AUC <- (K091_no_yes_day1$AUC)/K091_no_no_day1$meanAUC[1]
K147_no_yes_day1$relative_AUC <- (K147_no_yes_day1$AUC)/K147_no_no_day1$meanAUC[1]
K153_no_yes_day1$relative_AUC <- (K153_no_yes_day1$AUC)/K153_no_no_day1$meanAUC[1]
K163_no_yes_day1$relative_AUC <- (K163_no_yes_day1$AUC)/K163_no_no_day1$meanAUC[1]
K209_no_yes_day1$relative_AUC <- (K209_no_yes_day1$AUC)/K209_no_no_day1$meanAUC[1]

H53_no_yes_day3$relative_AUC <- (H53_no_yes_day3$AUC)/H53_no_no_day1$meanAUC[1]
K091_no_yes_day3$relative_AUC <- (K091_no_yes_day3$AUC)/K091_no_no_day1$meanAUC[1]
K147_no_yes_day3$relative_AUC <- (K147_no_yes_day3$AUC)/K147_no_no_day1$meanAUC[1]
K153_no_yes_day3$relative_AUC <- (K153_no_yes_day3$AUC)/K153_no_no_day1$meanAUC[1]
K163_no_yes_day3$relative_AUC <- (K163_no_yes_day3$AUC)/K163_no_no_day1$meanAUC[1]
K209_no_yes_day3$relative_AUC <- (K209_no_yes_day3$AUC)/K209_no_no_day1$meanAUC[1]

H53_no_yes_day5$relative_AUC <- (H53_no_yes_day5$AUC)/H53_no_no_day1$meanAUC[1]
K091_no_yes_day5$relative_AUC <- (K091_no_yes_day5$AUC)/K091_no_no_day1$meanAUC[1]
K147_no_yes_day5$relative_AUC <- (K147_no_yes_day5$AUC)/K147_no_no_day1$meanAUC[1]
K153_no_yes_day5$relative_AUC <- (K153_no_yes_day5$AUC)/K153_no_no_day1$meanAUC[1]
K163_no_yes_day5$relative_AUC <- (K163_no_yes_day5$AUC)/K163_no_no_day1$meanAUC[1]
K209_no_yes_day5$relative_AUC <- (K209_no_yes_day5$AUC)/K209_no_no_day1$meanAUC[1]

H53_no_yes_day15$relative_AUC <- (H53_no_yes_day15$AUC)/H53_no_no_day1$meanAUC[1]
K091_no_yes_day15$relative_AUC <- (K091_no_yes_day15$AUC)/K091_no_no_day1$meanAUC[1]
K147_no_yes_day15$relative_AUC <- (K147_no_yes_day15$AUC)/K147_no_no_day1$meanAUC[1]
K153_no_yes_day15$relative_AUC <- (K153_no_yes_day15$AUC)/K153_no_no_day1$meanAUC[1]
K163_no_yes_day15$relative_AUC <- (K163_no_yes_day15$AUC)/K163_no_no_day1$meanAUC[1]
K209_no_yes_day15$relative_AUC <- (K209_no_yes_day15$AUC)/K209_no_no_day1$meanAUC[1]

H53_yes_yes_day1$relative_AUC <- (H53_yes_yes_day1$AUC)/H53_no_no_day1$meanAUC
K091_yes_yes_day1$relative_AUC <- (K091_yes_yes_day1$AUC)/K091_no_no_day1$meanAUC[1]
K147_yes_yes_day1$relative_AUC <- (K147_yes_yes_day1$AUC)/K147_no_no_day1$meanAUC
K153_yes_yes_day1$relative_AUC <- (K153_yes_yes_day1$AUC)/K153_no_no_day1$meanAUC
K163_yes_yes_day1$relative_AUC <- (K163_yes_yes_day1$AUC)/K163_no_no_day1$meanAUC[1]
K209_yes_yes_day1$relative_AUC <- (K209_yes_yes_day1$AUC)/K209_no_no_day1$meanAUC[1]

H53_yes_yes_day3$relative_AUC <- (H53_yes_yes_day3$AUC)/H53_no_no_day1$meanAUC
K091_yes_yes_day3$relative_AUC <- (K091_yes_yes_day3$AUC)/K091_no_no_day1$meanAUC[1]
K147_yes_yes_day3$relative_AUC <- (K147_yes_yes_day3$AUC)/K147_no_no_day1$meanAUC
K153_yes_yes_day3$relative_AUC <- (K153_yes_yes_day3$AUC)/K153_no_no_day1$meanAUC
K163_yes_yes_day3$relative_AUC <- (K163_yes_yes_day3$AUC)/K163_no_no_day1$meanAUC[1]
K209_yes_yes_day3$relative_AUC <- (K209_yes_yes_day3$AUC)/K209_no_no_day1$meanAUC[1]

H53_yes_yes_day5$relative_AUC <- (H53_yes_yes_day5$AUC)/H53_no_no_day1$meanAUC
K091_yes_yes_day5$relative_AUC <- (K091_yes_yes_day5$AUC)/K091_no_no_day1$meanAUC[1]
K147_yes_yes_day5$relative_AUC <- (K147_yes_yes_day5$AUC)/K147_no_no_day1$meanAUC
K153_yes_yes_day5$relative_AUC <- (K153_yes_yes_day5$AUC)/K153_no_no_day1$meanAUC
K163_yes_yes_day5$relative_AUC <- (K163_yes_yes_day5$AUC)/K163_no_no_day1$meanAUC[1]
K209_yes_yes_day5$relative_AUC <- (K209_yes_yes_day5$AUC)/K209_no_no_day1$meanAUC[1]

H53_yes_yes_day15$relative_AUC <- (H53_yes_yes_day15$AUC)/H53_no_no_day1$meanAUC
K091_yes_yes_day15$relative_AUC <- (K091_yes_yes_day15$AUC)/K091_no_no_day1$meanAUC[1]
K147_yes_yes_day15$relative_AUC <- (K147_yes_yes_day15$AUC)/K147_no_no_day1$meanAUC
K153_yes_yes_day15$relative_AUC <- (K153_yes_yes_day15$AUC)/K153_no_no_day1$meanAUC
K163_yes_yes_day15$relative_AUC <- (K163_yes_yes_day15$AUC)/K163_no_no_day1$meanAUC[1]
K209_yes_yes_day15$relative_AUC <- (K209_yes_yes_day15$AUC)/K209_no_no_day1$meanAUC[1]

# For the rest of the samples

C021_no_no_day1$relative_AUC <- (C021_no_no_day1$AUC)/C021_no_no_day1$meanAUC
C324_no_no_day1$relative_AUC <- (C324_no_no_day1$AUC)/C324_no_no_day1$meanAUC
C232_no_no_day1$relative_AUC <- (C232_no_no_day1$AUC)/C232_no_no_day1$meanAUC
C286_no_no_day1$relative_AUC <- (C286_no_no_day1$AUC)/C286_no_no_day1$meanAUC
C309_no_no_day1$relative_AUC <- (C309_no_no_day1$AUC)/C309_no_no_day1$meanAUC
K112_no_no_day1$relative_AUC <- (K112_no_no_day1$AUC)/K112_no_no_day1$meanAUC
CF12_no_no_day1$relative_AUC <- (CF12_no_no_day1$AUC)/CF12_no_no_day1$meanAUC
CF13_no_no_day1$relative_AUC <- (CF13_no_no_day1$AUC)/CF13_no_no_day1$meanAUC

C021_no_no_day3$relative_AUC <- (C021_no_no_day3$AUC)/C021_no_no_day1$meanAUC
C324_no_no_day3$relative_AUC <- (C324_no_no_day3$AUC)/C324_no_no_day1$meanAUC
C232_no_no_day3$relative_AUC <- (C232_no_no_day3$AUC)/C232_no_no_day1$meanAUC
C286_no_no_day3$relative_AUC <- (C286_no_no_day3$AUC)/C286_no_no_day1$meanAUC
C309_no_no_day3$relative_AUC <- (C309_no_no_day3$AUC)/C309_no_no_day1$meanAUC
K112_no_no_day3$relative_AUC <- (K112_no_no_day3$AUC)/K112_no_no_day1$meanAUC
CF12_no_no_day3$relative_AUC <- (CF12_no_no_day3$AUC)/CF12_no_no_day1$meanAUC
CF13_no_no_day3$relative_AUC <- (CF13_no_no_day3$AUC)/CF13_no_no_day1$meanAUC

C021_no_no_day5$relative_AUC <- (C021_no_no_day5$AUC)/C021_no_no_day1$meanAUC
C324_no_no_day5$relative_AUC <- (C324_no_no_day5$AUC)/C324_no_no_day1$meanAUC
C232_no_no_day5$relative_AUC <- (C232_no_no_day5$AUC)/C232_no_no_day1$meanAUC
C286_no_no_day5$relative_AUC <- (C286_no_no_day5$AUC)/C286_no_no_day1$meanAUC
C309_no_no_day5$relative_AUC <- (C309_no_no_day5$AUC)/C309_no_no_day1$meanAUC
K112_no_no_day5$relative_AUC <- (K112_no_no_day5$AUC)/K112_no_no_day1$meanAUC
CF12_no_no_day5$relative_AUC <- (CF12_no_no_day5$AUC)/CF12_no_no_day1$meanAUC
CF13_no_no_day5$relative_AUC <- (CF13_no_no_day5$AUC)/CF13_no_no_day1$meanAUC

C021_no_no_day15$relative_AUC <- (C021_no_no_day15$AUC)/C021_no_no_day1$meanAUC
C324_no_no_day15$relative_AUC <- (C324_no_no_day15$AUC)/C324_no_no_day1$meanAUC
C232_no_no_day15$relative_AUC <- (C232_no_no_day15$AUC)/C232_no_no_day1$meanAUC
C286_no_no_day15$relative_AUC <- (C286_no_no_day15$AUC)/C286_no_no_day1$meanAUC
C309_no_no_day15$relative_AUC <- (C309_no_no_day15$AUC)/C309_no_no_day1$meanAUC
K112_no_no_day15$relative_AUC <- (K112_no_no_day15$AUC)/K112_no_no_day1$meanAUC
CF12_no_no_day15$relative_AUC <- (CF12_no_no_day15$AUC)/CF12_no_no_day1$meanAUC
CF13_no_no_day15$relative_AUC <- (CF13_no_no_day15$AUC)/CF13_no_no_day1$meanAUC

# For samples with pOXA-48 in presence/absence of antibiotic

C021_no_yes_day1$relative_AUC <- (C021_no_yes_day1$AUC)/C021_no_no_day1$meanAUC[1]
C324_no_yes_day1$relative_AUC <- (C324_no_yes_day1$AUC)/C324_no_no_day1$meanAUC
C232_no_yes_day1$relative_AUC <- (C232_no_yes_day1$AUC)/C232_no_no_day1$meanAUC
C286_no_yes_day1$relative_AUC <- (C286_no_yes_day1$AUC)/C286_no_no_day1$meanAUC
C309_no_yes_day1$relative_AUC <- (C309_no_yes_day1$AUC)/C309_no_no_day1$meanAUC[1]
K112_no_yes_day1$relative_AUC <- (K112_no_yes_day1$AUC)/K112_no_no_day1$meanAUC[1]
CF12_no_yes_day1$relative_AUC <- (CF12_no_yes_day1$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_no_yes_day1$relative_AUC <- (CF13_no_yes_day1$AUC)/CF13_no_no_day1$meanAUC[1]

C021_no_yes_day3$relative_AUC <- (C021_no_yes_day3$AUC)/C021_no_no_day1$meanAUC[1]
C324_no_yes_day3$relative_AUC <- (C324_no_yes_day3$AUC)/C324_no_no_day1$meanAUC
C232_no_yes_day3$relative_AUC <- (C232_no_yes_day3$AUC)/C232_no_no_day1$meanAUC
C286_no_yes_day3$relative_AUC <- (C286_no_yes_day3$AUC)/C286_no_no_day1$meanAUC
C309_no_yes_day3$relative_AUC <- (C309_no_yes_day3$AUC)/C309_no_no_day1$meanAUC[1]
K112_no_yes_day3$relative_AUC <- (K112_no_yes_day3$AUC)/K112_no_no_day1$meanAUC[1]
CF12_no_yes_day3$relative_AUC <- (CF12_no_yes_day3$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_no_yes_day3$relative_AUC <- (CF13_no_yes_day3$AUC)/CF13_no_no_day1$meanAUC[1]

C021_no_yes_day5$relative_AUC <- (C021_no_yes_day5$AUC)/C021_no_no_day1$meanAUC[1]
C324_no_yes_day5$relative_AUC <- (C324_no_yes_day5$AUC)/C324_no_no_day1$meanAUC
C232_no_yes_day5$relative_AUC <- (C232_no_yes_day5$AUC)/C232_no_no_day1$meanAUC
C286_no_yes_day5$relative_AUC <- (C286_no_yes_day5$AUC)/C286_no_no_day1$meanAUC
C309_no_yes_day5$relative_AUC <- (C309_no_yes_day5$AUC)/C309_no_no_day1$meanAUC[1]
K112_no_yes_day5$relative_AUC <- (K112_no_yes_day5$AUC)/K112_no_no_day1$meanAUC[1]
CF12_no_yes_day5$relative_AUC <- (CF12_no_yes_day5$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_no_yes_day5$relative_AUC <- (CF13_no_yes_day5$AUC)/CF13_no_no_day1$meanAUC[1]

C021_no_yes_day15$relative_AUC <- (C021_no_yes_day15$AUC)/C021_no_no_day1$meanAUC[1]
C324_no_yes_day15$relative_AUC <- (C324_no_yes_day15$AUC)/C324_no_no_day1$meanAUC
C232_no_yes_day15$relative_AUC <- (C232_no_yes_day15$AUC)/C232_no_no_day1$meanAUC
C286_no_yes_day15$relative_AUC <- (C286_no_yes_day15$AUC)/C286_no_no_day1$meanAUC
C309_no_yes_day15$relative_AUC <- (C309_no_yes_day15$AUC)/C309_no_no_day1$meanAUC[1]
K112_no_yes_day15$relative_AUC <- (K112_no_yes_day15$AUC)/K112_no_no_day1$meanAUC[1]
CF12_no_yes_day15$relative_AUC <- (CF12_no_yes_day15$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_no_yes_day15$relative_AUC <- (CF13_no_yes_day15$AUC)/CF13_no_no_day1$meanAUC[1]

C021_yes_yes_day1$relative_AUC <- (C021_yes_yes_day1$AUC)/C021_no_no_day1$meanAUC
C324_yes_yes_day1$relative_AUC <- (C324_yes_yes_day1$AUC)/C324_no_no_day1$meanAUC[1]
C232_yes_yes_day1$relative_AUC <- (C232_yes_yes_day1$AUC)/C232_no_no_day1$meanAUC
C286_yes_yes_day1$relative_AUC <- (C286_yes_yes_day1$AUC)/C286_no_no_day1$meanAUC[1]
C309_yes_yes_day1$relative_AUC <- (C309_yes_yes_day1$AUC)/C309_no_no_day1$meanAUC[1]
K112_yes_yes_day1$relative_AUC <- (K112_yes_yes_day1$AUC)/K112_no_no_day1$meanAUC[1]
CF12_yes_yes_day1$relative_AUC <- (CF12_yes_yes_day1$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_yes_yes_day1$relative_AUC <- (CF13_yes_yes_day1$AUC)/CF13_no_no_day1$meanAUC[1]

C021_yes_yes_day3$relative_AUC <- (C021_yes_yes_day3$AUC)/C021_no_no_day1$meanAUC
C324_yes_yes_day3$relative_AUC <- (C324_yes_yes_day3$AUC)/C324_no_no_day1$meanAUC[1]
C232_yes_yes_day3$relative_AUC <- (C232_yes_yes_day3$AUC)/C232_no_no_day1$meanAUC
C286_yes_yes_day3$relative_AUC <- (C286_yes_yes_day3$AUC)/C286_no_no_day1$meanAUC[1]
C309_yes_yes_day3$relative_AUC <- (C309_yes_yes_day3$AUC)/C309_no_no_day1$meanAUC[1]
K112_yes_yes_day3$relative_AUC <- (K112_yes_yes_day3$AUC)/K112_no_no_day1$meanAUC[1]
CF12_yes_yes_day3$relative_AUC <- (CF12_yes_yes_day3$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_yes_yes_day3$relative_AUC <- (CF13_yes_yes_day3$AUC)/CF13_no_no_day1$meanAUC[1]

C021_yes_yes_day5$relative_AUC <- (C021_yes_yes_day5$AUC)/C021_no_no_day1$meanAUC
C324_yes_yes_day5$relative_AUC <- (C324_yes_yes_day5$AUC)/C324_no_no_day1$meanAUC[1]
C232_yes_yes_day5$relative_AUC <- (C232_yes_yes_day5$AUC)/C232_no_no_day1$meanAUC
C286_yes_yes_day5$relative_AUC <- (C286_yes_yes_day5$AUC)/C286_no_no_day1$meanAUC[1]
C309_yes_yes_day5$relative_AUC <- (C309_yes_yes_day5$AUC)/C309_no_no_day1$meanAUC[1]
K112_yes_yes_day5$relative_AUC <- (K112_yes_yes_day5$AUC)/K112_no_no_day1$meanAUC[1]
CF12_yes_yes_day5$relative_AUC <- (CF12_yes_yes_day5$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_yes_yes_day5$relative_AUC <- (CF13_yes_yes_day5$AUC)/CF13_no_no_day1$meanAUC[1]

C021_yes_yes_day15$relative_AUC <- (C021_yes_yes_day15$AUC)/C021_no_no_day1$meanAUC
C324_yes_yes_day15$relative_AUC <- (C324_yes_yes_day15$AUC)/C324_no_no_day1$meanAUC[1]
C232_yes_yes_day15$relative_AUC <- (C232_yes_yes_day15$AUC)/C232_no_no_day1$meanAUC
C286_yes_yes_day15$relative_AUC <- (C286_yes_yes_day15$AUC)/C286_no_no_day1$meanAUC[1]
C309_yes_yes_day15$relative_AUC <- (C309_yes_yes_day15$AUC)/C309_no_no_day1$meanAUC[1]
K112_yes_yes_day15$relative_AUC <- (K112_yes_yes_day15$AUC)/K112_no_no_day1$meanAUC[1]
CF12_yes_yes_day15$relative_AUC <- (CF12_yes_yes_day15$AUC)/CF12_no_no_day1$meanAUC[1]
CF13_yes_yes_day15$relative_AUC <- (CF13_yes_yes_day15$AUC)/CF13_no_no_day1$meanAUC[1]

# Now, include everything in a table to furtherly represent it

df_list <- list(H53_no_no_day1, H53_no_no_day3, H53_no_no_day5, H53_no_no_day15,
                K091_no_no_day1, K091_no_no_day3, K091_no_no_day5, K091_no_no_day15,
                K147_no_no_day1, K147_no_no_day3, K147_no_no_day5, K147_no_no_day15,
                K153_no_no_day1, K153_no_no_day3, K153_no_no_day5, K153_no_no_day15,
                K163_no_no_day1, K163_no_no_day3, K163_no_no_day5, K163_no_no_day15,
                K209_no_no_day1, K209_no_no_day3, K209_no_no_day5, K209_no_no_day15,
                H53_no_yes_day1, H53_no_yes_day3, H53_no_yes_day5, H53_no_yes_day15,
                K091_no_yes_day1, K091_no_yes_day3, K091_no_yes_day5, K091_no_yes_day15,
                K147_no_yes_day1, K147_no_yes_day3, K147_no_yes_day5, K147_no_yes_day15,
                K153_no_yes_day1, K153_no_yes_day3, K153_no_yes_day5, K153_no_yes_day15,
                K163_no_yes_day1, K163_no_yes_day3, K163_no_yes_day5, K163_no_yes_day15,
                K209_no_yes_day1, K209_no_yes_day3, K209_no_yes_day5, K209_no_yes_day15,
                H53_yes_yes_day1, H53_yes_yes_day3, H53_yes_yes_day5, H53_yes_yes_day15,
                K091_yes_yes_day1, K091_yes_yes_day3, K091_yes_yes_day5, K091_yes_yes_day15,
                K147_yes_yes_day1, K147_yes_yes_day3, K147_yes_yes_day5, K147_yes_yes_day15,
                K153_yes_yes_day1, K153_yes_yes_day3, K153_yes_yes_day5, K153_yes_yes_day15,
                K163_yes_yes_day1, K163_yes_yes_day3, K163_yes_yes_day5, K163_yes_yes_day15,
                K209_yes_yes_day1, K209_yes_yes_day3, K209_yes_yes_day5, K209_yes_yes_day15,
                C021_no_no_day1, C021_no_no_day3, C021_no_no_day5, C021_no_no_day15,
                C324_no_no_day1, C324_no_no_day3, C324_no_no_day5, C324_no_no_day15,
                C232_no_no_day1, C232_no_no_day3, C232_no_no_day5, C232_no_no_day15,
                C286_no_no_day1, C286_no_no_day3, C286_no_no_day5, C286_no_no_day15,
                C309_no_no_day1, C309_no_no_day3, C309_no_no_day5, C309_no_no_day15,
                K112_no_no_day1, K112_no_no_day3, K112_no_no_day5, K112_no_no_day15,
                CF12_no_no_day1, CF12_no_no_day3, CF12_no_no_day5, CF12_no_no_day15,
                CF13_no_no_day1, CF13_no_no_day3, CF13_no_no_day5, CF13_no_no_day15,
                C021_no_yes_day1, C021_no_yes_day3, C021_no_yes_day5, C021_no_yes_day15,
                C324_no_yes_day1, C324_no_yes_day3, C324_no_yes_day5, C324_no_yes_day15,
                C232_no_yes_day1, C232_no_yes_day3, C232_no_yes_day5, C232_no_yes_day15,
                C286_no_yes_day1, C286_no_yes_day3, C286_no_yes_day5, C286_no_yes_day15,
                C309_no_yes_day1, C309_no_yes_day3, C309_no_yes_day5, C309_no_yes_day15,
                K112_no_yes_day1, K112_no_yes_day3, K112_no_yes_day5, K112_no_yes_day15,
                CF12_no_yes_day1, CF12_no_yes_day3, CF12_no_yes_day5, CF12_no_yes_day15,
                CF13_no_yes_day1, CF13_no_yes_day3, CF13_no_yes_day5, CF13_no_yes_day15,
                C021_yes_yes_day1, C021_yes_yes_day3, C021_yes_yes_day5, C021_yes_yes_day15,
                C324_yes_yes_day1, C324_yes_yes_day3, C324_yes_yes_day5, C324_yes_yes_day15,
                C232_yes_yes_day1, C232_yes_yes_day3, C232_yes_yes_day5, C232_yes_yes_day15,
                C286_yes_yes_day1, C286_yes_yes_day3, C286_yes_yes_day5, C286_yes_yes_day15,
                C309_yes_yes_day1, C309_yes_yes_day3, C309_yes_yes_day5, C309_yes_yes_day15,
                K112_yes_yes_day1, K112_yes_yes_day3, K112_yes_yes_day5, K112_yes_yes_day15,
                CF12_yes_yes_day1, CF12_yes_yes_day3, CF12_yes_yes_day5, CF12_yes_yes_day15,
                CF13_yes_yes_day1, CF13_yes_yes_day3, CF13_yes_yes_day5, CF13_yes_yes_day15)

whole_AUC_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Now, represent the relative values grouping by strain, depending on condition (WITH pOXA-48 and presence/absence of Ab) and grouped by day

plot_A <- whole_AUC_df %>%
  filter(Sample == "H53", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="H53 without Ab"); plot_A

plot_B <- whole_AUC_df %>%
  filter(Sample == "H53", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="H53 with Ab"); plot_B

plot_C <- whole_AUC_df %>%
  filter(Sample == "K091", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K091 without Ab"); plot_C

plot_D <- whole_AUC_df %>%
  filter(Sample == "K091", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K091 with Ab"); plot_D

plot_E <- whole_AUC_df %>%
  filter(Sample == "K147", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K147 without Ab"); plot_E

plot_F <- whole_AUC_df %>%
  filter(Sample == "K147", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K147 with Ab"); plot_F

plot_G <- whole_AUC_df %>%
  filter(Sample == "K153", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K153 without Ab"); plot_G

plot_H <- whole_AUC_df %>%
  filter(Sample == "K153", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K153 with Ab"); plot_H

plot_I <- whole_AUC_df %>%
  filter(Sample == "K163", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K163 without Ab"); plot_I

plot_J <- whole_AUC_df %>%
  filter(Sample == "K163", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K163 with Ab"); plot_J

plot_K <- whole_AUC_df %>%
  filter(Sample == "K209", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K209 without Ab"); plot_K

plot_L <- whole_AUC_df %>%
  filter(Sample == "K209", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="K209 with Ab"); plot_L

plot_M <- whole_AUC_df %>%
  filter(Sample == "C021", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C021 without Ab"); plot_M

plot_N <- whole_AUC_df %>%
  filter(Sample == "C021", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C021 with Ab"); plot_N

plot_O <- whole_AUC_df %>%
  filter(Sample == "C232", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C232 without Ab"); plot_O

plot_P <- whole_AUC_df %>%
  filter(Sample == "C232", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C232 with Ab"); plot_P

plot_Q <- whole_AUC_df %>%
  filter(Sample == "C286", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C286 without Ab"); plot_Q

plot_R <- whole_AUC_df %>%
  filter(Sample == "C286", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C286 with Ab"); plot_R

plot_S <- whole_AUC_df %>%
  filter(Sample == "C309", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C309 without Ab"); plot_S

plot_T <- whole_AUC_df %>%
  filter(Sample == "C309", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C309 with Ab"); plot_T

plot_U <- whole_AUC_df %>%
  filter(Sample == "C324", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C324 without Ab"); plot_U

plot_V <- whole_AUC_df %>%
  filter(Sample == "C324", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="C324 with Ab"); plot_V

plot_W <- whole_AUC_df %>%
  filter(Sample == "CF12", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="CF12 without Ab"); plot_W

plot_X <- whole_AUC_df %>%
  filter(Sample == "CF12", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="CF12 with Ab"); plot_X

plot_Y <- whole_AUC_df %>%
  filter(Sample == "CF13", Plasmid == "YES", Antibiotic == "NO") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="CF13 without Ab"); plot_Y

plot_Z <- whole_AUC_df %>%
  filter(Sample == "CF13", Plasmid == "YES", Antibiotic == "YES") %>%
  ggplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_boxplot(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_point(aes(y = relative_AUC, x = Day, color = Day, fill = Day, group = Day)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw() +
  labs(title="CF13 with Ab"); plot_Z

plot_A + plot_B + plot_C + plot_D + plot_E + plot_F + plot_G + plot_H + plot_I + plot_J + plot_K + plot_L +
  plot_M + plot_N + plot_O + plot_P + plot_Q + plot_R + plot_S + plot_T + plot_U + plot_V + plot_W + plot_X + plot_Y + plot_Z

plot_A + plot_B + plot_C + plot_D + plot_E + plot_F + plot_G + plot_H + plot_I + plot_J + plot_K + plot_L +
  plot_M + plot_N + plot_O + plot_P + plot_Q + plot_R + plot_S + plot_T + plot_U + plot_V + plot_W + plot_X + plot_Y + plot_Z

plot_W + plot_X + plot_Y + plot_Z + plot_A + plot_B + plot_C + plot_D + plot_E + plot_F + plot_G + plot_H + plot_I + plot_J + plot_K + plot_L

# Plot costs at day 1
#filtering outliers

whole_AUC_df_2 <- whole_AUC_df %>%
  filter(Sample != "C021", Sample != "C232", Sample != "C309", Sample != "C324", Sample != "C286", Sample != "K112")

whole_AUC_df_2 <- whole_AUC_df_2[-c(403,475),] # delete failed values

whole_AUC_df_2 %>%
  filter(Plasmid == "YES" & Antibiotic == "NO" & Day == 1) %>%
  ggplot(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_boxplot(aes(y = relative_AUC, x = Sample, color = Sample), lwd=0.9) +
  geom_point(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_hline(yintercept=1,linetype="longdash", color="darkgrey", size = 1.2) +
  theme_bw(base_size=20) + ylim(0.75, 1.3)

whole_AUC_df_2 %>%
  filter(Plasmid == "YES" & Antibiotic == "YES" & Day == 1) %>%
  ggplot(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_boxplot(aes(y = relative_AUC, x = Sample, color = Sample), lwd=0.9) +
  geom_point(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_hline(yintercept=1,linetype="longdash", color="darkgrey", size = 1.2) +
  theme_bw(base_size=20) + ylim(0.75, 1.3)

# Costs at day 3

whole_AUC_df_2 <- whole_AUC_df

whole_AUC_df_2 %>%
  filter(Plasmid == "YES" & Antibiotic == "NO" & Day == 3) %>%
  ggplot(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_boxplot(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_point(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw()

whole_AUC_df_2 %>%
  filter(Plasmid == "YES" & Antibiotic == "YES" & Day == 3) %>%
  ggplot(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_boxplot(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_point(aes(y = relative_AUC, x = Sample, color = Sample)) +
  geom_hline(yintercept=1, linetype="dotted", color="darkgrey") +
  theme_bw()