"""
TSS Anchoring Analysis
Maps the absolute positional distance of the structurally constrained 
20 bp and 24 bp motif pairs relative to the Transcription Start Site (TSS).
"""

import os
import re
import pandas as pd
import numpy as np

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa", "sorghum", "fonio"]
    
    # Target configurations
    targets = {
        "20bp_bZIP_Dof": r"(ACGT.{20}AAAG|AAAG.{20}ACGT|ACGT.{20}CTTT|CTTT.{20}ACGT)",
        "24bp_Dof_Dof": r"(AAAG.{24}AAAG|CTTT.{24}CTTT)" # Focus on the tandem functional orientations
    }
    
    positional_data = {"20bp_bZIP_Dof": [], "24bp_Dof_Dof": []}
    
    print("Mapping motif pair coordinates relative to the TSS...")
    
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
            seq_length = len(seq)
            
            for name, pattern in targets.items():
                for match in re.finditer(pattern, seq):
                    # Calculate distance from the 3' end (TSS) to the start of the motif complex
                    dist_to_tss = seq_length - match.start()
                    positional_data[name].append(dist_to_tss)

    print("\n=== TSS ANCHORING RESULTS ===\n")
    for name, distances in positional_data.items():
        if distances:
            median_dist = np.median(distances)
            mean_dist = np.mean(distances)
            print(f"[{name} Complex]")
            print(f"  Total mapped complexes: {len(distances)}")
            print(f"  Median distance to TSS: {median_dist:.1f} bp")
            print(f"  Mean distance to TSS:   {mean_dist:.1f} bp")
            
            # Check for proximal clustering (e.g., within 500bp)
            proximal_count = sum(1 for d in distances if d <= 500)
            proximal_percent = (proximal_count / len(distances)) * 100
            print(f"  Anchored within 500 bp of TSS: {proximal_percent:.1f}%\n")

if __name__ == "__main__":
    main()
