import os
import re
import gzip

crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]
base_dir = os.path.expanduser("~/promoter_project")
gene2go_path = os.path.join(base_dir, "gene2go.gz")

stress_go_terms = {
    "GO:0006950", "GO:0009628", "GO:0009414", 
    "GO:0009408", "GO:0009607"
}

print("Loading NCBI master GO database...")
stress_gene_ids = set()

# Parse the gene2go file to find all numerical GeneIDs associated with stress
with gzip.open(gene2go_path, 'rt') as f:
    next(f) # Skip header
    for line in f:
        columns = line.strip().split('\t')
        if len(columns) > 2:
            gene_id = columns[1]
            go_id = columns[2]
            if go_id in stress_go_terms:
                stress_gene_ids.add(gene_id)

print(f"Identified {len(stress_gene_ids)} total stress-associated GeneIDs globally.")
print("Mapping back to crop-specific locus tags...")

for crop in crops:
    gff3_path = os.path.join(base_dir, crop, "annotation.gff3")
    output_path = os.path.join(base_dir, crop, "stress_genes.txt")
    
    if not os.path.exists(gff3_path):
        continue
        
    crop_stress_genes = set()
    
    with open(gff3_path, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
                
            columns = line.strip().split('\t')
            if len(columns) < 9 or columns[2] != "gene":
                continue
                
            attributes = columns[8]
            
            # Extract the numerical GeneID (e.g., GeneID:839580)
            geneid_match = re.search(r'GeneID:(\d+)', attributes)
            if geneid_match:
                extracted_geneid = geneid_match.group(1)
                
                # If this GeneID is in our global stress list, grab the locus tag
                if extracted_geneid in stress_gene_ids:
                    # Look for locus_tag=AT1G01010 or ID=gene-AT1G01010
                    locus_match = re.search(r'locus_tag=([^;]+)', attributes)
                    if not locus_match:
                        locus_match = re.search(r'ID=gene-([^;]+)', attributes)
                        
                    if locus_match:
                        crop_stress_genes.add(locus_match.group(1))
                        
    with open(output_path, 'w') as out:
        for gene in sorted(crop_stress_genes):
            out.write(f"{gene}\n")
            
    print(f"Extracted {len(crop_stress_genes)} stress genes for {crop}.")

print("\nExtraction complete. Ready for statistical enrichment analysis.")
