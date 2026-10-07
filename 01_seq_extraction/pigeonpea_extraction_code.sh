mkdir -p ~/promoter_project/pigeon_pea
cd ~/promoter_project/pigeon_pea

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/340/665/GCF_000340665.1_C.cajan_V1.0/GCF_000340665.1_C.cajan_V1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/340/665/GCF_000340665.1_C.cajan_V1.0/GCF_000340665.1_C.cajan_V1.0_genomic.gff.gz

gunzip GCF_000340665.1_C.cajan_V1.0_genomic.fna.gz
gunzip GCF_000340665.1_C.cajan_V1.0_genomic.gff.gz

mv GCF_000340665.1_C.cajan_V1.0_genomic.fna genome.fasta
mv GCF_000340665.1_C.cajan_V1.0_genomic.gff annotation.gff3

samtools faidx genome.fasta
cut -f1,2 genome.fasta.fai > genome.chrom.sizes
awk -v OFS='\t' '$3=="gene" {print $1, $4-1, $5, $9, ".", $7}' annotation.gff3 > genes.bed
bedtools flank -i genes.bed -g genome.chrom.sizes -l 1000 -r 0 -s > promoters_coords.bed
bedtools getfasta -fi genome.fasta -bed promoters_coords.bed -s -nameOnly > true_promoters.fasta

ls -lh true_promoters.fasta
head -n 5 true_promoters.fasta