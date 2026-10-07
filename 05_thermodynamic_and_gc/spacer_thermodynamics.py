"""
Comprehensive Spacer Thermodynamics
Calculates localized GC content of the exact spacer sequences between 
all four motif orientations across 1 to 30 bp distances.
"""

import os
import re
import pandas as pd
from collections import defaultdict

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa", "sorghum", "fonio"]
    pairs = [("ACGT", "ACGT"), ("AAAG", "AAAG"), ("ACGT", "AAAG"), ("AAAG", "ACGT")]
    
    spacer_gc_data = defaultdict(list)
    
    print("Extracting localized spacer thermodynamics for all combinations (1-30 bp)...")
    
    for crop in crops:
        fasta_file = os.path.join(base_dir, crop, "true_promoters.fasta")
        if not os.path.exists(fasta_file):
            continue
            
        with open(fasta_file, 'r') as f:
            content = f.read().split('>')[1:]
            
        for entry in content:
            lines = entry.strip().split('\n')
            if not lines: continue
            seq = "".join(lines[1:]).upper()
            
            for m1, m2 in pairs:
                pair_name = f"{m1}_{m2}"
                for dist in range(1, 31):
                    # Using lookahead (?=...) to ensure overlapping pairs are caught
                    pattern = f"(?=({m1}(.{{{dist}}}){m2}))"
                    for match in re.finditer(pattern, seq):
                        spacer = match.group(2) 
                        gc = ((spacer.count('G') + spacer.count('C')) / len(spacer)) * 100
                        spacer_gc_data[f"{pair_name},{dist}"].append(gc)

    summary_stats = []
    for key, gc_list in spacer_gc_data.items():
        pair, dist = key.split(',')
        avg_gc = sum(gc_list) / len(gc_list)
        summary_stats.append({
            "Motif_Pair": pair,
            "Distance_bp": int(dist),
            "Spacer_Count": len(gc_list),
            "Avg_Spacer_GC": round(avg_gc, 2)
        })
        
    df = pd.DataFrame(summary_stats)
    
    if not df.empty:
        df = df.sort_values(by=["Motif_Pair", "Distance_bp"])
        
    output_csv = os.path.join(base_dir, "spacer_thermodynamics_summary.csv")
    df.to_csv(output_csv, index=False)
    
    print(f"\nAnalysis complete! Processed {sum(len(v) for v in spacer_gc_data.values())} total functional spacers.")
    print(f"Results saved to: {output_csv}")

if __name__ == "__main__":
    main()
