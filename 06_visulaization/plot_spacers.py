import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set style for publication-ready figures
plt.rcParams['figure.dpi'] = 300
sns.set_theme(style="whitegrid", context="paper")

base_dir = os.path.expanduser("~/promoter_project")
input_csv = os.path.join(base_dir, "rigorous_spatial_grammar.csv")
output_png = os.path.join(base_dir, "spatial_grammar_histograms.png")

print("Loading data...")
df = pd.read_csv(input_csv)

# Helper function to extract and flatten all distances from a specific column
def extract_distances(column_name):
    distances = []
    # Drop empty rows and filter out "NA"
    valid_rows = df[column_name].dropna()
    valid_rows = valid_rows[valid_rows != "NA"]
    
    for row in valid_rows:
        # Split by ';' and convert to integers
        distances.extend([int(x) for x in str(row).split(';')])
    return distances

# Define the four pairings to plot
pairings = {
    'ACGT followed by ACGT': 'ACGT_ACGT_Spacers',
    'ACGT followed by AAAG': 'ACGT_AAAG_Spacers',
    'AAAG followed by AAAG': 'AAAG_AAAG_Spacers',
    'AAAG followed by ACGT': 'AAAG_ACGT_Spacers'
}

# Create a 2x2 figure grid
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()
colors = sns.color_palette("deep", 4)

print("Plotting distributions...")
for idx, (title, col) in enumerate(pairings.items()):
    ax = axes[idx]
    dists = extract_distances(col)
    
    if len(dists) > 0:
        # Plot histogram with 1-bp bins up to 30bp, overlay with a density curve (KDE)
        sns.histplot(dists, bins=range(0, 32), kde=True, color=colors[idx], ax=ax, edgecolor='black', alpha=0.7)
        ax.set_title(f"{title}\n(Total pairs = {len(dists)})", fontsize=12, fontweight='bold')
        ax.set_ylabel("Frequency", fontsize=11)
        
        # Highlight potential helical peaks (e.g., 10bp, 20bp) if they exist
        ax.axvline(10, color='red', linestyle='--', alpha=0.5, label='1 Turn (~10bp)')
        ax.axvline(20, color='red', linestyle=':', alpha=0.5, label='2 Turns (~20bp)')
        if idx == 0:
            ax.legend(loc='upper right')
    else:
        ax.set_title(f"{title}\n(No pairs found)", fontsize=12, fontweight='bold')

    # X-axis labels only on the bottom row
    if idx >= 2:
        ax.set_xlabel("Spacer Distance (bp)", fontsize=11)
    
    ax.set_xlim(0, 31)

plt.tight_layout()
plt.savefig(output_png, bbox_inches='tight')
print(f"Success! High-resolution figure saved to: {output_png}")

