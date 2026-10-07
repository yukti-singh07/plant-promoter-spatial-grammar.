mkdir -p ~/promoter_project/arabidopsis
cd ~/promoter_project/arabidopsis

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip GCF_000001735.4_TAIR10.1_genomic.gff.gz

mv GCF_000001735.4_TAIR10.1_genomic.fna genome.fasta
mv GCF_000001735.4_TAIR10.1_genomic.gff annotation.gff3

samtools faidx genome.fasta
cut -f1,2 genome.fasta.fai > genome.chrom.sizes
awk -v OFS='\t' '$3=="gene" {print $1, $4-1, $5, $9, ".", $7}' annotation.gff3 > genes.bed
bedtools flank -i genes.bed -g genome.chrom.sizes -l 1000 -r 0 -s > promoters_coords.bed
bedtools getfasta -fi genome.fasta -bed promoters_coords.bed -s -nameOnly > true_promoters.fasta

ls -lh true_promoters.fasta
head -n 5 true_promoters.fasta