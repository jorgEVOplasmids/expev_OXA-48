# Analysis of the evolution of plasmid-mediated antibiotic resistance through WGS

## Introduction

In this repository you'll find the code developed for the analysis of WGS data resulting from the experimental evolution of clinical bacterial strains carrying pOXA-48, a carbapenem-resistance plasmid of great relevance.

All the scripts included in this project have been developed in order to automatize already developed programs, and to parse, cure and depict all the information obtained by multiple analysis on the evolved bacteria during the experiment (variant calling, analysis of growth curves, analysis of plasmid copy number, etc.).

This README.md file is a brief guide that explains the utility of each script developed during the project. The structure of the repository follows the same as the materials and methods section in the work for the sake of simplicity.

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

## From raw reads to variant calling

