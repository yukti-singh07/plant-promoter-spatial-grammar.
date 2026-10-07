"""
Promoter Sequence Environment Analyzer
Calculates GC content and extracts motif flanking sequences (5bp) 
from stress-responsive promoters across all crops.
"""

import os
import re
import pandas as pd

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]

    gc_results = []
    acgt_flanks = []
    aaag_flanks = []

    print("Analyzing sequence environments for stress promoters...")

    for crop in crops:
        stress_file = os.path.join(base_dir, crop, "stress_genes.txt")
        fasta_file = os.path.join(base_dir, crop, "true_promoters.fasta")
        
        if not os.path.exists(stress_file) or not os.path.exists(fasta_file):
            continue
            
        # 1. Load the stress genes for this crop
        stress_genes = set()
        with open(stress_file, 'r') as f:
            for line in f:
                if line.strip():
                    stress_genes.add(line.strip())
                    
        # 2. Parse the FASTA file manually (no Biopython required)
        sequences = {}
        with open(fasta_file, 'r') as f:
            header = ""
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        sequences[header] = "".join(seq)
                    header = line[1:] 
                    seq = []
                else:
                    seq.append(line.upper())
            if header:
                sequences[header] = "".join(seq)
                
        # 3. Process the stress-responsive promoters
        for header, seq in sequences.items():
            is_stress = False
            matched_gene = ""
            
            # Match the FASTA header to our stress gene list
            for gene in stress_genes:
                if gene in header:
                    is_stress = True
                    matched_gene = gene
                    break
            
            if is_stress:
                # A. Calculate GC Content
                g_count = seq.count('G')
                c_count = seq.count('C')
                length = len(seq)
                gc_percent = ((g_count + c_count) / length) * 100 if length > 0 else 0
                
                gc_results.append({
                    "Species": crop,
                    "Gene_ID": matched_gene,
                    "Promoter_Length": length,
                    "GC_Percent": round(gc_percent, 2)
                })
                
                # B. Extract Flanking sequences (5bp on each side)
                # Find all ACGT occurrences
                for match in re.finditer(r'ACGT', seq):
                    start, end = match.start(), match.end()
                    if start >= 5 and end + 5 <= len(seq):
                        flank_seq = seq[start-5:end+5]
                        acgt_flanks.append(f">{crop}_{matched_gene}_ACGT\n{flank_seq}")
                        
                # Find all AAAG occurrences
                for match in re.finditer(r'AAAG', seq):
                    start, end = match.start(), match.end()
                    if start >= 5 and end + 5 <= len(seq):
                        flank_seq = seq[start-5:end+5]
                        aaag_flanks.append(f">{crop}_{matched_gene}_AAAG\n{flank_seq}")

    # 4. Save Outputs
    gc_out = os.path.join(base_dir, "stress_gc_content.csv")
    pd.DataFrame(gc_results).to_csv(gc_out, index=False)
    
    acgt_out = os.path.join(base_dir, "acgt_flanks_for_weblogo.fasta")
    with open(acgt_out, "w") as f:
        f.write("\n".join(acgt_flanks) + "\n")
        
    aaag_out = os.path.join(base_dir, "aaag_flanks_for_weblogo.fasta")
    with open(aaag_out, "w") as f:
        f.write("\n".join(aaag_flanks) + "\n")
        
    print(f"\nSaved GC content to: {gc_out}")
    print(f"Saved ACGT sequences to: {acgt_out}")
    print(f"Saved AAAG sequences to: {aaag_out}")

if __name__ == "__main__":
    main()
