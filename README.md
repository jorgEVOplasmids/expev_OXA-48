# Analysis of the evolution of plasmid-mediated antibiotic resistance through WGS

## Introduction

In this repository you'll find the code developed for the analysis of WGS data resulting from the experimental evolution of clinical bacterial strains carrying pOXA-48, a carbapenem-resistance plasmid of great relevance.

All the scripts included in this project have been developed in order to automatize the execuution of already developed programs, and to parse, cure and depict all the information obtained by multiple analysis on the evolved bacteria during the experiment (variant calling, analysis of growth curves, analysis of plasmid copy number, etc.).

This README.md file is a brief guide that explains the utility of each script developed during the project. The structure of the repository follows the same as the Materials and Methods section in the work for the sake of simplicity.

Additionally to these Python and R scripts, other software have been necessary to thoroughly analyze this data:

* FastQC (https://github.com/s-andrews/FastQC)
* MultiQC (https://github.com/ewels/MultiQC)
* BBMap Clumpify (https://github.com/BioInfoTools/BBMap)
* Trim Galore (https://github.com/FelixKrueger/TrimGalore)
* SPAdes (https://github.com/ablab/spades)
* QUAST (https://github.com/ablab/quast)
* BWA-MEM (https://github.com/lh3/bwa)
* SAMtools (https://github.com/samtools/samtools)
* GNU datamash (https://www.gnu.org/software/datamash/)
* Qualimap (http://qualimap.conesalab.org/)
* Breseq (https://github.com/barricklab/breseq)
* Snippy (https://github.com/tseemann/snippy)
* iJump (https://github.com/sleyn/ijump)
* panISa (https://github.com/bvalot/panISa)
* ISfinder (https://isfinder.biotoul.fr/)
* Anvi'o (https://anvio.org/)
* Prokka (https://github.com/tseemann/prokka)

![fig_workflow (1)-1](https://user-images.githubusercontent.com/105753869/169649481-1a4c46b9-41c7-4f42-a443-ce5c48b9f90b.png)

*Summary scheme of the bioinformatic workflow carried out during this project. Input/intermediate data is shown in blue; output data is shown in green. Most scripts included in this repository were useful to process all these outputs given by already developed programs mentioned above.*

## Analysis of growth curves

Growth curves data was analyzed using R. Multiple factors were studied. Firstly, fitness of each plasmid-carrying sample (with and without antibiotic) was compared to its respective plasmid-free sample fitness. This relative fitness gives information about **plasmid cost**. Moreover, growth rate at days 3,5 and 15 of each plasmid-carrying sample was compared to its respective plasmid-free to check possible **cost compensation**. Finally, bacterial growth at the end and beginning of the EE was compared to see if there had been an **improvement in fitness during the evolution**. All these comparisons are performed in the same R script (**analyses_curves_stats_and_plots.R**), as well as the code for plotting all the results.

File:

https://github.com/jorgEVOplasmids/expev_pOXA-48/blob/main/curves/analyses_curves_stats_and_plots.R

## From raw reads to variant calling

In order to automatize the execution of all the programs necessary for carrying out multiple steps of the analysis (Processing of raw reads + *De novo* assembly + Mapping against reference + Variant calling) the script **breseq_pipeline.sh** was developed. It wraps all the commands necessary for their execution from the UNIX terminal with the specific parameters indicated during the Materials and Methods section.

File:

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/rawreads_to_vcall/breseq_pipeline.sh

## Variant calling

### Parsing breseq output

As breseq outputs an HTML file per sample analyzed, the Python script **breseq_index_parser.py** was developed to transform the HTML output files from all the replicates and samples of each strain and merge them in a single XLSX so that their mutations could be easily compared.

File: 

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/variant_calling/breseq_index_parser.py

### Analyzing mutation spectrum

To study the mutation spectrum in the different EE groups, the XLSX files obtained in the previous step were processed with another Python script (**mut_spectrum.py**). It basically checks the results from breseq and summarizes them depending on the nature of the mutations and their type (synonymous, non-synonymous or intergenic; SNP/Indel or New Junction).

File:

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/variant_calling/mut_spectrum.py

### Summary and visualization of variant calling results

For examining parallel evolution events, breseq results were displayed as heatmaps showing the frequency of mutations affecting each gene in each strain replicate depending on its EE group and mutation type (SNPs/Indels or NJ evidences). For this, firstly, the frequency of all mutations affecting a gene were summarized for each strain using **mut_screening.py**. However, breseq results indicated the mutations affecting genes with the PGAP annotation. Many of these do not count with a gene name and are annotated with a generic identifier (pgap_ + number). Thus, for analyzing the variant calling results in an easier way, those genes annotated with the PGAP generic identifier were renamed with the PGAP description (which gives information about the gene function and sometimes includes the actual gene name) with the script **rename_genes.py**.

Then, for displaying the information of multiple samples together, specifically, of *K. pneumoniae* samples, the data was merged and displayed with **merge_klebsiella.py**.

Files:

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/variant_calling/mut_screening.py

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/variant_calling/rename_genes.py

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/variant_calling/merge_klebsiella.py

## Analysis of IS rearrangements

### Summary of IS movements

To study the influence of IS movements in each EE group, a Python script which detects the variant calling events annotated as IS was written. Besides, the calculous of the percentage of IS movements due to the pOXA-48 was included in the code. Finally, statistical analysis comparing the number of IS rearrangements in presence and in absence of pOXA-48 were performed in R.

Files:

### Identification of pOXA-48 IS targets

For identifying the genes targeted by pOXA-48 IS elements, **is_targets.py** was written. It parses the XSLX files generated by breseq_index_parser.py and filters those genes which appear mutated by a NJ including a pOXA-48 IS.

File:

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/is_analysis/is_targets.py

## Analysis of PCN

### Statistical analysis of PCN and plotting

For the statistical analyses carried out on PCN data (multiple comparisons and LME), and representing the PCN of all the plasmids contained in each strain replicate, the script **pcn_plots_stats.R** was developed. It includes all the statistics and code for plotting the boxplots shown in the project.

File:



## Visualization and functional analysis with Anvi'o

### Running anvi'o populations workflow

All the anvi'o programs included in the populations workflow indicated by its developers were included in a script for automatizing the execution of all the necessary steps (https://merenlab.org/2015/07/20/analyzing-variability/).

File: 

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/anvio/anvio-populations.sh
