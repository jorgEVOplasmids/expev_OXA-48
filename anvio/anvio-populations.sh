#!/bin/bash

### Basic microbial population genetics pipeline with Anvi'o ###

# Allow activation of conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Iterate through strains folders
for folder in reads_trimmed/*/
do
	strain=$( basename $folder )
	
	# Annotate samples with prokka
	conda activate prokka
	prokka --proteins ref_genomes/$strain*.gbk --outdir anvio/$strain/functions/ --prefix $strain.annot ref_genomes/$strain.reformat.fasta --cpus 5
	conda deactivate
	conda activate anvio-7.1
	
	# Get necessary txt files
	python ~/Documents/TFM/scripts/gff_parser.py anvio/$strain/functions/$strain.annot.gff --gene-calls anvio/$strain/functions/gene_calls.txt --annotation anvio/$strain/functions/gene_annot.txt
	
	# Create an Anvi'o contigs db from closed bacteria reference sequence and import annotations
	anvi-gen-contigs-database -T 5 -L 5000 -f ref_genomes/$strain.reformat.fasta -o anvio/$strain/complete_contigs.db --external-gene-calls anvio/$strain/functions/gene_calls.txt -n 'Contigs db from closed '$strain' reference'
	
	# Import functions
	anvi-import-functions -c anvio/$strain/complete_contigs.db -i anvio/$strain/functions/gene_annot.txt
	
	# Create profile db for each sample
	for bamfile in mapping/$strain/*.sorted.bam
	do
		samplename=$( basename $bamfile )
		samplename=${samplename::-11}
		anvi-profile -T 5 -i $bamfile -c anvio/$strain/complete_contigs.db -o anvio/$strain/profiles/$samplename-profiledb/ -S $samplename
	done
	
	# Merge all profile databases into one
	profiles=$( echo anvio/$strain/profiles/*/PROFILE.db )
	anvi-merge  $profiles -o anvio/$strain/profiles/$strain-merged -c anvio/$strain/complete_contigs.db
done

conda deactivate

# For visualization execute next line with the strain desired
# anvi-interactive -p anvio/$strain/profiles/$strain-merged/PROFILE.db -c anvio/$strain/complete_contigs.db

