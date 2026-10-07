"""
Generates a publication-ready grouped bar chart for macro promoter motif percentages.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    input_csv = os.path.join(base_dir, "macro_all_promoters_results.csv")
    output_fig = os.path.join(base_dir, "macro_motif_distribution.png")

    if not os.path.exists(input_csv):
        print("Results CSV not found. Run the macro analysis script first.")
        return

    df = pd.read_csv(input_csv)

    # Set overall aesthetic style
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 6))

    # Create grouped bar plot
    ax = sns.barplot(
        data=df,
        x="Species",
        y="Percentage_%",
        hue="Motif_Pair",
        palette="muted"
    )

    # Formatting the plot
    plt.title("Global Macro-Level Distribution of Cis-Regulatory Pairs (1–20 bp Window)", fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Species", fontsize=12, fontweight='bold')
    plt.ylabel("Promoters Containing Motif Pair (%)", fontsize=12, fontweight='bold')
    plt.xticks(rotation=30, ha='right', fontsize=10)
    plt.legend(title="Motif Pair", title_fontsize='11', fontsize='10', bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.tight_layout()

    # Save as high-res PNG and PDF for publication
    plt.savefig(output_fig, dpi=300)
    plt.savefig(output_fig.replace(".png", ".pdf"))
    print(f"\nFigure successfully saved to {output_fig} (and .pdf)")

if __name__ == "__main__":
    main()
