import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    csv_dir = os.path.expanduser("~/promoter_project/go_results")
    if not os.path.exists(csv_dir):
        print(f"Directory {csv_dir} not found.")
        return

    all_data = []
    for filepath in glob.glob(os.path.join(csv_dir, "*.csv")):
        filename = os.path.basename(filepath).replace(".csv", "")
        try:
            df = pd.read_csv(filepath)
            if "source" in df.columns:
                df = df[df["source"] == "GO:BP"]
            df = df.sort_values(by="adjusted_p_value").head(5)
            for _, row in df.iterrows():
                all_data.append({
                    "Architecture": filename,
                    "Pathway": row["term_name"],
                    "Gene_Count": row["intersection_size"],
                    "P_adj": row["adjusted_p_value"],
                    "MinusLog10P": row["negative_log10_of_adjusted_p_value"]
                })
        except Exception as e:
            print(f"Error processing {filepath}: {e}")

    if not all_data:
        print("No data found.")
        return

    combined_df = pd.DataFrame(all_data)
    combined_df["Spacer"] = combined_df["Architecture"].str.extract(r"(\d+)").astype(int)
    combined_df = combined_df.sort_values(by=["Spacer", "Architecture"])

    plt.figure(figsize=(12, 8))
    sns.set_theme(style="whitegrid")
    
    sns.scatterplot(
        data=combined_df,
        x="Architecture",
        y="Pathway",
        size="Gene_Count",
        hue="MinusLog10P",
        palette="viridis",
        sizes=(100, 800),
        alpha=0.8,
        edgecolor="black"
    )
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0., title="-log10(Padj) / Count")
    plt.xticks(rotation=45, ha="right", fontsize=12, fontweight="bold")
    plt.yticks(fontsize=11)
    plt.xlabel("Spatial Architecture", fontsize=14, fontweight="bold", labelpad=15)
    plt.ylabel("Biological Process (GO:BP)", fontsize=14, fontweight="bold", labelpad=15)
    plt.title("Functional Divergence Across Promoter Architectures", fontsize=16, fontweight="bold", pad=20)
    
    plt.tight_layout()
    output_img = os.path.join(csv_dir, "Figure_5_Functional_Bubble_Plot.png")
    plt.savefig(output_img, dpi=300, bbox_inches="tight")
    print(f"Success! Plot saved to {output_img}")

if __name__ == "__main__":
    main()
