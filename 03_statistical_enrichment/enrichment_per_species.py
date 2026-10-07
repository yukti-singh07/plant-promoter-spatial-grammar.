"""
Species-Specific Spatial Grammar Enrichment Analysis
Calculates Fisher's Exact statistics for motif spacer distances 
individually for each crop using the rigorous AAAG core motif.
"""

import os
import re
import pandas as pd
import scipy.stats as stats

def extract_locus(attr_string):
    match = re.search(r'locus_tag=([^;]+)', str(attr_string))
    if not match:
        match = re.search(r'ID=gene-([^;]+)', str(attr_string))
    return match.group(1) if match else str(attr_string)

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    input_csv = os.path.join(base_dir, "rigorous_spatial_grammar.csv")
    output_csv = os.path.join(base_dir, "enrichment_per_species.csv")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]

    print("Loading master spatial grammar dataset...")
    df = pd.read_csv(input_csv)
    df['Clean_Gene_ID'] = df['Gene_ID'].apply(extract_locus)

    # REVERTED DICTIONARY: Using the strict 4-letter AAAG motif
    pairings = {
        'ACGT_ACGT': 'ACGT_ACGT_Spacers',
        'ACGT_AAAG': 'ACGT_AAAG_Spacers',
        'AAAG_AAAG': 'AAAG_AAAG_Spacers',
        'AAAG_ACGT': 'AAAG_ACGT_Spacers'
    }

    results = []

    for crop in crops:
        print(f"Processing {crop}...")
        
        # Isolate the whole genome background for ONLY this crop
        crop_df = df[df['Species'] == crop].copy()
        if crop_df.empty:
            continue
            
        # Load the stress genes for ONLY this crop
        stress_file = os.path.join(base_dir, crop, "stress_genes.txt")
        stress_set = set()
        if os.path.exists(stress_file):
            with open(stress_file, 'r') as f:
                for line in f:
                    gene_id = line.strip()
                    if gene_id:
                        stress_set.add(gene_id)
        
        crop_df['Is_Stress'] = crop_df['Clean_Gene_ID'].isin(stress_set)
        
        if crop_df['Is_Stress'].sum() == 0:
            continue

        for pair_name, col_name in pairings.items():
            if col_name not in crop_df.columns:
                print(f"  Warning: Column {col_name} not found in master CSV. Skipping.")
                continue
                
            valid_rows = crop_df.dropna(subset=[col_name])
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

                results.append({
                    "Species": crop,
                    "Motif_Pair": pair_name,
                    "Distance_bp": dist,
                    "Stress_Genes_With_Dist": stress_with,
                    "Background_Genes_With_Dist": bg_with,
                    "Odds_Ratio": round(odds_ratio, 4) if odds_ratio != float('inf') else "Inf",
                    "Fisher_P_Value": fisher_p,
                    "Significant_Enrichment": "YES" if fisher_p < 0.05 and odds_ratio > 1 else "NO"
                })

    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df.to_csv(output_csv, index=False)
        print(f"\nSpecies-specific enrichment complete. Saved to {output_csv}")
        
        significant = results_df[results_df["Significant_Enrichment"] == "YES"]
        print(f"Found {len(significant)} significant spacer enrichments across individual species.")
    else:
        print("\nNo results generated. Ensure your master CSV contains the spacer columns.")

if __name__ == "__main__":
    main()
