"""
Spatial Grammar Enrichment Analysis
Calculates Fisher's Exact and Chi-Square statistics for motif spacer 
distances (0-30 bp) in stress-responsive vs. background promoters.
"""

import os
import re
import pandas as pd
import scipy.stats as stats

def extract_locus(attr_string):
    """Extracts the clean gene name from a messy GFF3 attribute string."""
    match = re.search(r'locus_tag=([^;]+)', str(attr_string))
    if not match:
        match = re.search(r'ID=gene-([^;]+)', str(attr_string))
    return match.group(1) if match else str(attr_string)

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    input_csv = os.path.join(base_dir, "rigorous_spatial_grammar.csv")
    output_csv = os.path.join(base_dir, "enrichment_results.csv")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]

    print("Loading master spatial grammar dataset...")
    df = pd.read_csv(input_csv)

    print("Cleaning GFF3 attribute strings...")
    df['Clean_Gene_ID'] = df['Gene_ID'].apply(extract_locus)

    # 1. Load Stress Genes and create a universal identifier
    stress_set = set()
    for crop in crops:
        stress_file = os.path.join(base_dir, crop, "stress_genes.txt")
        if os.path.exists(stress_file):
            with open(stress_file, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    if gene_id:
                        stress_set.add(f"{crop}_{gene_id}")
    
    # 2. Tag each row in the dataframe as Stress (True) or Background (False)
    df['Universal_ID'] = df['Species'] + '_' + df['Clean_Gene_ID']
    df['Is_Stress'] = df['Universal_ID'].isin(stress_set)
    
    print(f"Identified {df['Is_Stress'].sum()} total stress-responsive promoters matching the CSV.")

    # 3. Setup statistical tracking
    pairings = {
        'ACGT_ACGT': 'ACGT_ACGT_Spacers',
        'ACGT_AAAG': 'ACGT_AAAG_Spacers',
        'AAAG_AAAG': 'AAAG_AAAG_Spacers',
        'AAAG_ACGT': 'AAAG_ACGT_Spacers'
    }
    
    results = []

    # 4. Run tests for every distance from 0 to 30 bp
    print("Calculating Chi-Square and Fisher's Exact tests...")
    for pair_name, col_name in pairings.items():
        valid_rows = df.dropna(subset=[col_name])
        valid_rows = valid_rows[valid_rows[col_name] != "NA"].copy()
        
        if valid_rows.empty:
            continue

        valid_rows['Parsed_Spacers'] = valid_rows[col_name].apply(lambda x: [int(d) for d in str(x).split(';')])
        
        for dist in range(31):
            valid_rows['Has_Distance'] = valid_rows['Parsed_Spacers'].apply(lambda spacers: dist in spacers)
            
            stress_with = len(valid_rows[(valid_rows['Is_Stress'] == True) & (valid_rows['Has_Distance'] == True)])
            stress_without = len(valid_rows[(valid_rows['Is_Stress'] == True) & (valid_rows['Has_Distance'] == False)])
            bg_with = len(valid_rows[(valid_rows['Is_Stress'] == False) & (valid_rows['Has_Distance'] == True)])
            bg_without = len(valid_rows[(valid_rows['Is_Stress'] == False) & (valid_rows['Has_Distance'] == False)])
            
            table = [[stress_with, stress_without], [bg_with, bg_without]]
            
            if stress_with + bg_with == 0:
                continue
                
            odds_ratio, fisher_p = stats.fisher_exact(table, alternative='two-sided')
            
            try:
                chi2, chi2_p, dof, expected = stats.chi2_contingency(table)
            except ValueError:
                chi2, chi2_p = "NA", "NA"

            results.append({
                "Motif_Pair": pair_name,
                "Distance_bp": dist,
                "Stress_Genes_With_Dist": stress_with,
                "Background_Genes_With_Dist": bg_with,
                "Odds_Ratio": round(odds_ratio, 4) if odds_ratio != float('inf') else "Inf",
                "Fisher_P_Value": fisher_p,
                "Chi2_P_Value": chi2_p,
                "Significant_Enrichment": "YES" if fisher_p < 0.05 and odds_ratio > 1 else "NO"
            })

    # 5. Save Results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv, index=False)
    print(f"\nEnrichment analysis complete. Statistics saved to {output_csv}")
    
    significant = results_df[results_df["Significant_Enrichment"] == "YES"]
    print(f"\nFound {len(significant)} statistically significant spacer enrichments across all pairs.")

if __name__ == "__main__":
    main()