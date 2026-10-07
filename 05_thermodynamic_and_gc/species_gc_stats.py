"""
Advanced GC Content Statistical Analysis - Table Generator (Species-Wise)
Compares stress vs. background GC content individually for each crop 
plus pooled clade and overall summaries using exact gene ID extraction.
"""

import os
import re
import pandas as pd
from scipy.stats import mannwhitneyu

def extract_locus(header):
    """Extracts the clean gene ID from a GFF3-derived FASTA header.
    Handles both locus_tag= and ID=gene- formats used across NCBI RefSeq species."""
    match = re.search(r'locus_tag=([^;]+)', header)
    if not match:
        match = re.search(r'ID=gene-([^;]+)', header)
    return match.group(1) if match else header.split()[0]

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    
    dicots = ["chickpea", "arabidopsis", "pigeon_pea", "cassava", "quinoa"]
    monocots = ["sorghum"]  # Fonio excluded due to missing annotation data
    all_crops = dicots + monocots
    
    data = []
    stats_summary = []

    print("Extracting genome-wide GC content with exact gene ID matching...")

    for crop in all_crops:
        clade = "Dicot" if crop in dicots else "Monocot"
        stress_file = os.path.join(base_dir, crop, "stress_genes.txt")
        fasta_file = os.path.join(base_dir, crop, "true_promoters.fasta")
        
        if not os.path.exists(fasta_file):
            print(f"Warning: Missing promoter file for {crop}, skipping.")
            continue
            
        stress_genes = set()
        if os.path.exists(stress_file):
            with open(stress_file, 'r') as f:
                for line in f:
                    if line.strip():
                        stress_genes.add(line.strip())
                        
        with open(fasta_file, 'r') as f:
            header = ""
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        seq_str = "".join(seq)
                        length = len(seq_str)
                        if length > 0:
                            gc = ((seq_str.count('G') + seq_str.count('C')) / length) * 100
                            gene_id = extract_locus(header)
                            is_stress = gene_id in stress_genes
                            data.append({"Species": crop, "Clade": clade, "Type": "Stress" if is_stress else "Background", "GC_Percent": gc})
                    header = line[1:]
                    seq = []
                else:
                    seq.append(line.upper())
            
            if header:
                seq_str = "".join(seq)
                length = len(seq_str)
                if length > 0:
                    gc = ((seq_str.count('G') + seq_str.count('C')) / length) * 100
                    gene_id = extract_locus(header)
                    is_stress = gene_id in stress_genes
                    data.append({"Species": crop, "Clade": clade, "Type": "Stress" if is_stress else "Background", "GC_Percent": gc})

    df = pd.DataFrame(data)
    
    def run_statistics(group_df, label):
        stress = group_df[group_df['Type'] == 'Stress']['GC_Percent']
        background = group_df[group_df['Type'] == 'Background']['GC_Percent']
        
        if len(stress) == 0 or len(background) == 0:
            return
            
        stat, p_value = mannwhitneyu(stress, background, alternative='two-sided')
        
        stats_summary.append({
            "Comparison_Group": label,
            "Stress_Count": len(stress),
            "Stress_Avg_GC": round(stress.mean(), 2),
            "Background_Count": len(background),
            "Background_Avg_GC": round(background.mean(), 2),
            "P_Value": f"{p_value:.2e}",
            "Significant": "YES" if p_value < 0.05 else "NO"
        })

    # 1. Run individual species comparisons
    for crop in all_crops:
        crop_df = df[df['Species'] == crop]
        run_statistics(crop_df, f"Species: {crop.capitalize()}")

    # 2. Run pooled and overall comparisons[cite: 1, 2]
    run_statistics(df[df['Clade'] == 'Dicot'], "POOLED: DICOTS ONLY")
    run_statistics(df[df['Clade'] == 'Monocot'], "POOLED: MONOCOTS ONLY")
    run_statistics(df, "OVERALL (All 6 Species)")
    
    # Save the raw data
    output_csv = os.path.join(base_dir, "comprehensive_gc_stats.csv")
    df.to_csv(output_csv, index=False)
    
    # Save the expanded summary table
    summary_csv = os.path.join(base_dir, "gc_statistics_summary.csv")
    pd.DataFrame(stats_summary).to_csv(summary_csv, index=False)
    
    print(f"\nSpecies-wise analysis complete! Summary table saved to: {summary_csv}")

if __name__ == "__main__":
    main()
