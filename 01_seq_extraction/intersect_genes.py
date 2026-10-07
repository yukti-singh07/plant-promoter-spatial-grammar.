import os

def main():
    base_dir = os.path.expanduser("~/promoter_project")
    crops = ["chickpea", "arabidopsis", "pigeon_pea", "sorghum", "cassava", "quinoa", "fonio"]
    
    motifs = [
        "1bp_bZIP_bZIP", "17bp_bZIP_bZIP", "26bp_bZIP_bZIP", "30bp_bZIP_bZIP", 
        "1bp_Dof_Dof", "8bp_Dof_Dof", "20bp_bZIP_Dof", "24bp_Dof_Dof"
    ]

    print("Loading stress-responsive baseline...")
    stress_master = set()
    
    for crop in crops:
        stress_file = os.path.join(base_dir, crop, "stress_genes.txt")
        if os.path.exists(stress_file):
            with open(stress_file, 'r') as f:
                genes = {line.strip() for line in f if line.strip()}
                stress_master.update(genes)
                
    print(f"Loaded {len(stress_master)} total stress-responsive gene IDs.")
    print("\nIntersecting with spatial architectures...")

    for motif in motifs:
        motif_file = os.path.join(base_dir, f"{motif}_GeneIDs.txt")
        if not os.path.exists(motif_file):
            continue
            
        with open(motif_file, 'r') as f:
            motif_genes = {line.strip() for line in f if line.strip()}
            
        active_stress_motif_genes = motif_genes.intersection(stress_master)
        
        output_file = os.path.join(base_dir, f"GO_READY_{motif}.txt")
        with open(output_file, 'w') as out:
            for gene in sorted(active_stress_motif_genes):
                out.write(f"{gene}\n")
                
        print(f"{motif}: {len(active_stress_motif_genes)} strictly filtered genes saved for GO analysis.")

if __name__ == "__main__":
    main()
