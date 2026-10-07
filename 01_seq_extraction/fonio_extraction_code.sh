mkdir -p ~/promoter_project/fonio
cd ~/promoter_project/fonio

# Download the Fonio CM05836 assembly and annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/859/565/GCA_902859565.1_Fonio_CM05836/GCA_902859565.1_Fonio_CM05836_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/859/565/GCA_902859565.1_Fonio_CM05836/GCA_902859565.1_Fonio_CM05836_genomic.gff.gz

# Decompress the files
gunzip GCA_902859565.1_Fonio_CM05836_genomic.fna.gz
gunzip GCA_902859565.1_Fonio_CM05836_genomic.gff.gz

# Rename files to standardize your workflow
mv GCA_902859565.1_Fonio_CM05836_genomic.fna genome.fasta
mv GCA_902859565.1_Fonio_CM05836_genomic.gff annotation.gff3

# Run the identical extraction pipeline
samtools faidx genome.fasta
cut -f1,2 genome.fasta.fai > genome.chrom.sizes
awk -v OFS='\t' '$3=="gene" {print $1, $4-1, $5, $9, ".", $7}' annotation.gff3 > genes.bed
bedtools flank -i genes.bed -g genome.chrom.sizes -l 1000 -r 0 -s > promoters_coords.bed
bedtools getfasta -fi genome.fasta -bed promoters_coords.bed -s -nameOnly > true_promoters.fasta

# Verify the final extraction
ls -lh true_promoters.fasta
head -n 5 true_promoters.fasta