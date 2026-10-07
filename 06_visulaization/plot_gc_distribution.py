"""
Publication Figure Generation: GC Content Distribution Violin Plot
Visualizes the spread of GC content between Stress and Background promoters
across the 6 individual species.
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    input_csv = os.path.join(base_dir, "comprehensive_gc_stats.csv")
    output_path = os.path.join(base_dir, "Figure_GC_Distribution.pdf")

    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found. Run your GC stats script first.")
        return

    df = pd.read_csv(input_csv)

    # Set up the high-resolution plot style
    fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
    
    # Create a clean seaborn boxplot or violin plot comparing Stress vs Background per species
    sns.boxplot(
        data=df, 
        x="Species", 
        y="GC_Percent", 
        hue="Type", 
        palette="Set2", 
        linewidth=1.2,
        fliersize=1,
        ax=ax
    )

    # Formatting for publication
    plt.title('Promoter GC Content Distribution: Stress-Responsive vs. Background Genes', fontsize=13, fontweight='bold', pad=15)
    plt.xlabel('Species', fontsize=11, fontweight='bold')
    plt.ylabel('Promoter GC Content (%)', fontsize=11, fontweight='bold')
    plt.xticks(rotation=20, fontsize=10)
    plt.legend(title='Promoter Type', title_fontsize='10', fontsize='9')
    
    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.savefig(output_path.replace(".pdf", ".png"), dpi=300, bbox_inches='tight')
    
    print(f"GC Distribution figure successfully saved to: {output_path}")

if __name__ == "__main__":
    main()
