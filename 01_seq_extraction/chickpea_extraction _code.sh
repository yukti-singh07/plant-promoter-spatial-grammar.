#!/bin/bash
# Chickpea (Cicer arietinum cv. CDC Frontier) Promoter Extraction Pipeline
# Assembly: GCF_000331145.1 (ASM33114v1)

# 1. Create and navigate to the working directory
mkdir -p ~/promoter_project/chickpea
cd ~/promoter_project/chickpea

# 2. Download the RefSeq genome sequence (FASTA) and annotation (GFF3) from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/331/145/GCF_000331145.1_ASM33114v1/GCF_000331145.1_ASM33114v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/331/145/GCF_000331145.1_ASM33114v1/GCF_000331145.1_ASM33114v1_genomic.gff.gz

# 3. Decompress the files
gunzip GCF_000331145.1_ASM33114v1_genomic.fna.gz
gunzip GCF_000331145.1_ASM33114v1_genomic.gff.gz

# 4. Rename files to standardize the workflow
mv GCF_000331145.1_ASM33114v1_genomic.fna genome.fasta
mv GCF_000331145.1_ASM33114v1_genomic.gff annotation.gff3

# 5. Index the genome and create a chromosome sizes file
samtools faidx genome.fasta
cut -f1,2 genome.fasta.fai > genome.chrom.sizes

# 6. Isolate 'gene' features from the GFF3 and format them into a 0-based BED file
awk -v OFS='\t' '$3=="gene" {print $1, $4-1, $5, $9, ".", $7}' annotation.gff3 > genes.bed

# 7. Calculate the exact -1000 to 0 bp upstream coordinates (accounting for strand direction)
bedtools flank -i genes.bed -g genome.chrom.sizes -l 1000 -r 0 -s > promoters_coords.bed

# 8. Extract the actual DNA sequences into the final FASTA file
bedtools getfasta -fi genome.fasta -bed promoters_coords.bed -s -nameOnly > true_promoters.fasta

# 9. Verify the output file
ls -lh true_promoters.fasta
head -n 5 true_promoters.fasta