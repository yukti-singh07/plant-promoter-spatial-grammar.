import os
import re

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]
    
    targets = {
        "1bp_bZIP_bZIP": r"(ACGT.{1}ACGT)",
        "17bp_bZIP_bZIP": r"(ACGT.{17}ACGT)",
        "26bp_bZIP_bZIP": r"(ACGT.{26}ACGT)",
        "30bp_bZIP_bZIP": r"(ACGT.{30}ACGT)",
        "1bp_Dof_Dof": r"(AAAG.{1}AAAG|CTTT.{1}CTTT)",
        "8bp_Dof_Dof": r"(AAAG.{8}AAAG|CTTT.{8}CTTT)",
        "20bp_bZIP_Dof": r"(ACGT.{20}AAAG|AAAG.{20}ACGT|ACGT.{20}CTTT|CTTT.{20}ACGT)",
        "24bp_Dof_Dof": r"(AAAG.{24}AAAG|CTTT.{24}CTTT)"
    }
    
    gene_lists = {name: set() for name in targets.keys()}
    
    for target_crop in crops:
        fasta_file = os.path.join(base_dir, target_crop, "true_promoters.fasta")
        
        if not os.path.exists(fasta_file):
            continue
            
        with open(fasta_file, 'r') as f:
            content = f.read().split('>')[1:]
            
        for entry in content:
            lines = entry.strip().split('\n')
            if not lines: continue
            
            header = lines[0]
            clean_id = None
            
            # 1. Catch Arabidopsis AGI codes
            agi_match = re.search(r'(AT[1-5CM]G\d{5})', header, re.IGNORECASE)
            # 2. Catch the GFF-style semicolon tags from your headers
            id_match = re.search(r'ID=gene-([^;]+)', header)
            gene_match = re.search(r'gene=([^;]+)', header)
            locus_match = re.search(r'locus_tag=([^;]+)', header)
            
            if agi_match:
                clean_id = agi_match.group(1).upper()
            elif id_match:
                clean_id = id_match.group(1)
            elif gene_match:
                clean_id = gene_match.group(1)
            elif locus_match:
                clean_id = locus_match.group(1)
            else:
                clean_id = header.split()[0]
                
            seq = "".join(lines[1:]).upper()
            
            for name, pattern in targets.items():
                if re.search(pattern, seq):
                    gene_lists[name].add(clean_id)

    for name, ids in gene_lists.items():
        output_file = os.path.join(base_dir, f"{name}_GeneIDs.txt")
        with open(output_file, 'w') as out_f:
            for gene_id in sorted(ids):
                out_f.write(f"{gene_id}\n")
                
        print(f"Extracted IDs for {name}")

if __name__ == "__main__":
    main()
