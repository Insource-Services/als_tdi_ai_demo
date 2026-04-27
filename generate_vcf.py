import pandas as pd
import numpy as np

def generate_synthetic_vcf():
    # 1. Load the dataset
    df = pd.read_csv('synthetic_proact_dataset.csv')
    n_patients = len(df)
    
    np.random.seed(42)
    
    # 2. Define variants
    variants = [
        {"id": "C9orf72_exp", "chrom": "chr9", "pos": 27573482, "ref": "C", "alt": "<EXP>"},
        {"id": "SOD1_A4V", "chrom": "chr21", "pos": 33039605, "ref": "C", "alt": "T"},
        {"id": "SOD1_D90A", "chrom": "chr21", "pos": 33031935, "ref": "A", "alt": "G"},
        {"id": "rs12486783", "chrom": "chr5", "pos": 151000000, "ref": "A", "alt": "G"}, # TNIP1
        {"id": "rs11674437", "chrom": "chr2", "pos": 100000000, "ref": "C", "alt": "T"}, # C2orf69
        {"id": "rs11736730", "chrom": "chr7", "pos": 40000000, "ref": "G", "alt": "A"},  # YKT6
        {"id": "rs35714649", "chrom": "chr17", "pos": 26000000, "ref": "T", "alt": "C"}  # SARM1
    ]
    
    # Initialize genotypes dict: key=variant_id, value=list of genotypes '0/0', '0/1', '1/1'
    genotypes = {v["id"]: ['0/0'] * n_patients for v in variants}
    
    # Identify patient indices
    bulbar_indices = df[df['Site_of_Onset'] == 'Bulbar'].index.tolist()
    limb_indices = df[df['Site_of_Onset'] == 'Limb'].index.tolist()
    
    # --- Constraint 1: Mendelian causal variants ---
    # C9orf72: n=80 bulbar onset
    c9_selected = np.random.choice(bulbar_indices, size=80, replace=False)
    for idx in c9_selected:
        genotypes["C9orf72_exp"][idx] = '0/1'
        
    # SOD1: n=40 limb onset (20 A4V, 20 D90A)
    sod1_selected = np.random.choice(limb_indices, size=40, replace=False)
    sod1_a4v = sod1_selected[:20]
    sod1_d90a = sod1_selected[20:]
    
    for idx in sod1_a4v:
        genotypes["SOD1_A4V"][idx] = '0/1'
    for idx in sod1_d90a:
        genotypes["SOD1_D90A"][idx] = '0/1'
        
    mendelian_indices = set(c9_selected).union(set(sod1_selected))
    sporadic_indices = list(set(range(n_patients)) - mendelian_indices)
    
    # --- Constraint 2: Non-mendelian (polygenic) risk alleles ---
    # YKT6: high risk rare variant assigned to 2% of sporadic cohort
    ykt6_n = int(len(sporadic_indices) * 0.02)
    ykt6_selected = np.random.choice(sporadic_indices, size=ykt6_n, replace=False)
    for idx in ykt6_selected:
        genotypes["rs11736730"][idx] = '0/1'
        
    # SARM1: correlated with rapid progression
    # Sort sporadic cohort by ALSFRS_R_Slope (ascending, most negative = fastest)
    sporadic_df = df.iloc[sporadic_indices].sort_values(by='ALSFRS_R_Slope')
    sorted_sporadic_idx = sporadic_df.index.tolist()
    
    # Assign higher probability to top 25% fastest progressors
    n_fast = int(len(sorted_sporadic_idx) * 0.25)
    fast_progressors = sorted_sporadic_idx[:n_fast]
    slow_progressors = sorted_sporadic_idx[n_fast:]
    
    for idx in fast_progressors:
        if np.random.rand() < 0.35: # 35% chance for fast progressors
            genotypes["rs35714649"][idx] = '0/1'
    for idx in slow_progressors:
        if np.random.rand() < 0.05: # 5% chance for slow progressors
            genotypes["rs35714649"][idx] = '0/1'
            
    # TNIP1 (rs12486783) MAF=0.38
    # C2orf69 (rs11674437) MAF=0.35
    # Assign across the entire cohort based on Hardy-Weinberg
    def assign_hw(maf, size):
        p_00 = (1 - maf)**2
        p_01 = 2 * maf * (1 - maf)
        p_11 = maf**2
        return np.random.choice(['0/0', '0/1', '1/1'], size=size, p=[p_00, p_01, p_11])
        
    tnip1_gts = assign_hw(0.38, n_patients)
    c2orf69_gts = assign_hw(0.35, n_patients)
    
    for i in range(n_patients):
        genotypes["rs12486783"][i] = tnip1_gts[i]
        genotypes["rs11674437"][i] = c2orf69_gts[i]
        
    # 3. Write to VCF file
    vcf_filename = "synthetic_cohort.vcf"
    with open(vcf_filename, 'w') as f:
        # Headers
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=SyntheticGenerator\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Name">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        # Column headers
        columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        patient_ids = df['Patient_ID'].tolist()
        f.write("\t".join(columns + patient_ids) + "\n")
        
        # Variants
        for v in variants:
            vid = v["id"]
            gene = ""
            if "C9orf72" in vid: gene = "C9orf72"
            elif "SOD1" in vid: gene = "SOD1"
            elif vid == "rs12486783": gene = "TNIP1"
            elif vid == "rs11674437": gene = "C2orf69"
            elif vid == "rs11736730": gene = "YKT6"
            elif vid == "rs35714649": gene = "SARM1"
            
            # Calculate simple AF for the INFO field
            gts = genotypes[vid]
            alt_count = sum(1 for gt in gts if '1' in gt) # naive count
            af = alt_count / (2 * n_patients)
            
            info = f"AF={af:.3f};GENE={gene}"
            
            row = [
                v["chrom"],
                str(v["pos"]),
                vid,
                v["ref"],
                v["alt"],
                ".",
                "PASS",
                info,
                "GT"
            ]
            row.extend(gts)
            f.write("\t".join(row) + "\n")
            
    print(f"Successfully generated {vcf_filename} for {n_patients} patients.")
    print("Mendelian and Polygenic constraints applied.")

if __name__ == "__main__":
    generate_synthetic_vcf()
