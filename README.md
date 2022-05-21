# Analysis of the evolution of plasmid-mediated antibiotic resistance through WGS

## Introduction

In this repository you'll find the code developed for the analysis of WGS data resulting from the experimental evolution of clinical bacterial strains carrying pOXA-48, a carbapenem-resistance plasmid of great relevance.

All the scripts included in this project have been developed in order to automatize already developed programs, and to parse, cure and depict all the information obtained by multiple analysis on the evolved bacteria during the experiment (variant calling, analysis of growth curves, analysis of plasmid copy number, etc.).

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

[fig_workflow_TFM.pdf](https://github.com/jorgEVOplasmids/expev_OXA-48/files/8747478/fig_workflow_TFM.pdf)

*Summary scheme of the bioinformatic workflow carried out during this project. Input/intermediate data is shown in blue; output data is shown in green. Most scripts included in this repository were useful to process all these outputs given by already developed programs mentioned above.*

## From raw reads to variant calling

In order to automatize the execution of all the programs necessary for carrying out multiple steps of the analysis (Processing of raw reads + *De novo* assembly + Mapping against reference + Variant calling) the script **breseq_pipeline.sh** was developed. It wraps all the commands necessary for their execution from the UNIX terminal with the specific parameters indicated during the Materials and Methods section.

## Variant calling

### Parsing breseq output

As breseq outputs an HTML file per sample analyzed, the Python script **breseq_index_parser.py** was developed to transform the HTML output files from all the replicates and samples of each strain and merge them in a single XLSX so that their mutations could be easily compared.

File: 

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/variant_calling/breseq_index_parser.py

### Analyzing mutation spectrum

To study the mutation spectrum in the different EE groups, the XLSX files obtained in the previous step were processed with another Python script (**mut_spectrum.py**). It basically checks the results from breseq and summarizes them depending on the nature of the mutations.

### Summary and visualization of variant calling results

For examining parallel evolution events, breseq results were displayed as heatmaps showing the frequency of mutations affecting each gene in each strain replicate depending on its EE group and mutation type (SNPs/Indels or NJ evidences). For this, firstly, the frequency of all mutations affecting a gene was summarized for each strain using **mut_screening.py**. Then, for displaying the information of multiple samples together, specifically, of *K. pneumoniae* samples, the data was merged and displayed with **merge_klebsiella.py**.

## Analysis of IS rearrangements

### Summary of IS movements

## Visualization and functional analysis with Anvi'o

### Running anvi'o populations workflow

All the anvi'o programs included in the populations workflow indicated by its developers were included in a script for automatizing the execution of all the necessary steps (https://merenlab.org/tutorials/363infant-gut/#chapter-vi-microbial-population-genetics).

File: 

https://github.com/jorgEVOplasmids/expev_OXA-48/blob/main/anvio/anvio-populations.sh
