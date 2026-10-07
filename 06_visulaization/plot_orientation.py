"""
Publication Figure Generation: Strand Orientation
Generates a grouped bar chart for the 20 bp and 24 bp complex orientations.
"""

import matplotlib.pyplot as plt
import numpy as np
import os

def main():
    # Set up the data
    labels = ['+/+ (Head-to-Tail)', '+/- (Head-to-Head)', '-/+ (Tail-to-Tail)', '-/- (Tail-to-Head)']
    
    # Exact counts from your statistical output
    dof_dof_24bp = [14263, 9966, 9299, 12756]
    bzip_dof_20bp = [2417, 2115, 2148, 2282]
    
    x = np.arange(len(labels))
    width = 0.35
    
    # Create the figure with high resolution
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Plot the bars
    rects1 = ax.bar(x - width/2, dof_dof_24bp, width, label='24 bp Dof-Dof (Flexible Hinge)', color='#2c7bb6')
    rects2 = ax.bar(x + width/2, bzip_dof_20bp, width, label='20 bp bZIP-Dof (Rigid Tether)', color='#d7191c')
    
    # Add text, labels, and custom x-axis tick labels
    ax.set_ylabel('Total Genome-wide Occurrences', fontsize=12, fontweight='bold')
    ax.set_title('Directional Docking Constraints of Spatially Enriched Motif Pairs', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.legend(fontsize=11)
    
    # Add statistical significance markers
    ax.text(0 - width/2, 14500, '***', ha='center', va='bottom', fontsize=14, fontweight='bold')
    ax.text(3 - width/2, 13000, '***', ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    # Clean up layout and save
    plt.tight_layout()
    output_path = os.path.expanduser("~/promoter_project/Figure_Orientation.pdf")
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    
    print(f"Figure successfully generated and saved to: {output_path}")

if __name__ == "__main__":
    main()
