"""
Motif Spatial Grammar Visualization Pipeline
Generates high-resolution spacer distance distributions (0-30 bp) for ACGT and AAAG pairings.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set style parameters for publication-ready figures
plt.rcParams['figure.dpi'] = 300
sns.set_theme(style="whitegrid", context="paper")

def main():
    # Define file paths (Users should update base_dir if running outside the standard project structure)
    base_dir = os.path.expanduser("~/promoter_project")
    input_csv = os.path.join(base_dir, "rigorous_spatial_grammar.csv")
    output_png = os.path.join(base_dir, "spatial_grammar_histograms.png")

    print(f"Loading spatial grammar data from {input_csv}...")
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print("Error: Input CSV not found. Please run analyze_motifs.py first.")
        return

    def extract_distances(column_name):
        """Extracts and flattens semicolon-separated distances from a specified column."""
        distances = []
        valid_rows = df[column_name].dropna()
        valid_rows = valid_rows[valid_rows != "NA"]
        
        for row in valid_rows:
            distances.extend([int(x) for x in str(row).split(';')])
        return distances

    # Map the target spatial orientations to their respective dataframe columns
    pairings = {
        'ACGT followed by ACGT': 'ACGT_ACGT_Spacers',
        'ACGT followed by AAAG': 'ACGT_AAAG_Spacers',
        'AAAG followed by AAAG': 'AAAG_AAAG_Spacers',
        'AAAG followed by ACGT': 'AAAG_ACGT_Spacers'
    }

    # Initialize a 2x2 plotting grid
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
    axes = axes.flatten()
    colors = sns.color_palette("deep", 4)

    print("Generating distribution plots...")
    for idx, (title, col) in enumerate(pairings.items()):
        ax = axes[idx]
        dists = extract_distances(col)
        
        if dists:
            # Plot the distance distribution with a 1-bp bin resolution and KDE curve
            sns.histplot(dists, bins=range(0, 32), kde=True, color=colors[idx], 
                         ax=ax, edgecolor='black', alpha=0.7)
            ax.set_title(f"{title}\n(Total pairs = {len(dists)})", fontsize=12, fontweight='bold')
            ax.set_ylabel("Frequency", fontsize=11)
            
            # Annotate potential DNA helical constraints (1 turn ~10bp, 2 turns ~20bp)
            ax.axvline(10, color='red', linestyle='--', alpha=0.5, label='1 Turn (~10bp)')
            ax.axvline(20, color='red', linestyle=':', alpha=0.5, label='2 Turns (~20bp)')
            
            if idx == 0:
                ax.legend(loc='upper right')
        else:
            ax.set_title(f"{title}\n(No pairs found)", fontsize=12, fontweight='bold')

        # Format X-axis labels exclusively on the bottom row to reduce clutter
        if idx >= 2:
            ax.set_xlabel("Spacer Distance (bp)", fontsize=11)
        
        # Enforce strict 0-30bp visualization window
        ax.set_xlim(0, 31)

    plt.tight_layout()
    plt.savefig(output_png, bbox_inches='tight')
    print(f"High-resolution figure successfully saved to: {output_png}")

if __name__ == "__main__":
    main()