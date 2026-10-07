"""
Publication Figure Generation: Evolutionary Enrichment Heatmap
Automatically calculates and visualizes the conservation of 20 bp and 24 bp 
spatial complexes across the 6 selected crop species from rigorous spatial data.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    input_csv = os.path.join(base_dir, "rigorous_spatial_grammar.csv")
    output_path = os.path.join(base_dir, "Figure_Enrichment_Heatmap.pdf")

    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found. Run your spatial grammar script first.")
        return

    df_data = pd.read_csv(input_csv)

    species_order = ["arabidopsis", "chickpea", "pigeon_pea", "cassava", "quinoa", "sorghum"]
    species_labels = [
        "Arabidopsis\n(Dicot)", "Chickpea\n(Dicot)", "Pigeon pea\n(Dicot)", 
        "Cassava\n(Dicot)", "Quinoa\n(Dicot)", "Sorghum\n(Monocot)"
    ]

    matrix_data = []
    complexes = ["24 bp Dof-Dof", "20 bp bZIP-Dof"]

    dof_dof_values = []
    bzip_dof_values = []

    for sp in species_order:
        sp_df = df_data[df_data["Species"] == sp]
        total_genes = len(sp_df)
        
        if total_genes == 0:
            dof_dof_values.append(0.0)
            bzip_dof_values.append(0.0)
            continue

        dof_dof_count = 0
        bzip_dof_count = 0

        for _, row in sp_df.iterrows():
            spacers_dof = str(row["AAAG_AAAG_Spacers"]).split(";")
            if any(s.strip() == "24" for s in spacers_dof):
                dof_dof_count += 1

            spacers_bzip1 = str(row["ACGT_AAAG_Spacers"]).split(";")
            spacers_bzip2 = str(row["AAAG_ACGT_Spacers"]).split(";")
            if any(s.strip() == "20" for s in spacers_bzip1) or any(s.strip() == "20" for s in spacers_bzip2):
                bzip_dof_count += 1

        dof_dof_values.append(round((dof_dof_count / total_genes) * 100, 2))
        bzip_dof_values.append(round((bzip_dof_count / total_genes) * 100, 2))

    matrix = np.array([dof_dof_values, bzip_dof_values])
    heat_df = pd.DataFrame(matrix, index=complexes, columns=species_labels)

    fig, ax = plt.subplots(figsize=(9, 4), dpi=300)

    sns.heatmap(
        heat_df, 
        annot=True,          
        fmt=".2f",           
        cmap="YlOrRd",       
        linewidths=1,        
        cbar_kws={'label': 'Promoters with Exact Spatial Lock (%)'} 
    )

    ax.axvline(5, color='black', linewidth=2, linestyle='--')

    plt.title('Evolutionary Conservation of Spatially Constrained Promoter Complexes', 
              fontsize=13, fontweight='bold', pad=15)
    plt.yticks(rotation=0, fontweight='bold')
    plt.xticks(fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.savefig(output_path.replace(".pdf", ".png"), dpi=300, bbox_inches='tight')
    
    print(f"\nHeatmap successfully generated and saved to: {output_path}")

if __name__ == "__main__":
    main()
