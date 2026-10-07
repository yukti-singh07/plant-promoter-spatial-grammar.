"""
Macro View: Percentage of ALL promoters containing motif pairs (1-20 bp window)
"""

import os
import re
from Bio import SeqIO
import pandas as pd

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]
    
    motifs = {
        'ACGT_ACGT': r'ACGT.{1,20}ACGT',
        'ACGT_AAAG': r'ACGT.{1,20}AAAG',
        'AAAG_AAAG': r'AAAG.{1,20}AAAG',
        'AAAG_ACGT': r'AAAG.{1,20}ACGT'
    }

    results = []

    for crop in crops:
        promoter_fasta = os.path.join(base_dir, crop, "true_promoters.fasta")
        if not os.path.exists(promoter_fasta):
            print(f"Skipping {crop}: promoters file not found.")
            continue
            
        promoters = list(SeqIO.parse(promoter_fasta, "fasta"))
        total_promoters = len(promoters)
        
        if total_promoters == 0:
            continue

        print(f"Processing {crop} ({total_promoters} total promoters)...")

        for pair_name, regex in motifs.items():
            count = 0
            for record in promoters:
                seq = str(record.seq).upper()
                if re.search(regex, seq):
                    count += 1
            
            percentage = (count / total_promoters) * 100
            results.append({
                "Species": crop,
                "Motif_Pair": pair_name.replace("_", "-"),
                "Total_Promoters": total_promoters,
                "Promoters_With_Pair": count,
                "Percentage_%": round(percentage, 2)
            })

    df = pd.DataFrame(results)
    output_csv = os.path.join(base_dir, "macro_all_promoters_results.csv")
    df.to_csv(output_csv, index=False)
    print(f"\nSaved clean macro promoter results to {output_csv}")
    print(df.to_string())

if __name__ == "__main__":
    main()
