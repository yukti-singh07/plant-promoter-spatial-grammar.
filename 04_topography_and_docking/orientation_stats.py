"""
Chi-Square Goodness-of-Fit Test for Strand Orientation
Validates whether the observed binding directions (+/+, +/-, -/+, -/-) 
deviate significantly from a random, uniform distribution.
"""

import scipy.stats as stats

def main():
    # Observed counts from the strand_orientation.py output
    # Order: [+/+ (Head-to-Tail), +/- (Head-to-Head), -/+ (Tail-to-Tail), -/- (Tail-to-Head)]
    
    dof_dof_24bp = [14263, 9966, 9299, 12756]
    bzip_dof_20bp = [2417, 2115, 2148, 2282] 
    
    print("=== CHI-SQUARE ORIENTATION STATISTICS ===\n")
    
    # 24 bp Analysis
    chi2_24, p_24 = stats.chisquare(f_obs=dof_dof_24bp)
    print("24 bp Homotypic (Dof-Dof) Complex:")
    print(f"Observed Counts: {dof_dof_24bp}")
    print(f"Chi-square statistic: {chi2_24:.2f}")
    print(f"p-value: {p_24:.2e}")
    if p_24 < 0.05:
        print("Verdict: SIGNIFICANT directional constraint (Tandem binding required).\n")
        
    # 20 bp Analysis
    chi2_20, p_20 = stats.chisquare(f_obs=bzip_dof_20bp)
    print("20 bp Heterotypic (bZIP-Dof) Complex:")
    print(f"Observed Counts: {bzip_dof_20bp}")
    print(f"Chi-square statistic: {chi2_20:.2f}")
    print(f"p-value: {p_20:.2e}")
    if chi2_20 < chi2_24:
        print("Verdict: HIGHLY PERMISSIVE relative to the 24 bp complex (Bidirectional binding supported).")

if __name__ == "__main__":
    main()