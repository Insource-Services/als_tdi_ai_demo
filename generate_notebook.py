import json

def create_notebook():
    notebook = {
        "cells": [],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "codemirror_mode": {"name": "ipython", "version": 3},
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython3",
                "version": "3.8.0"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 4
    }

    def add_markdown(text):
        notebook["cells"].append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [line + "\n" for line in text.split("\n")]
        })

    def add_code(code):
        notebook["cells"].append({
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [line + "\n" for line in code.split("\n")]
        })

    # Setup cell
    add_markdown("# ALS Genotype-Phenotype Analysis\nThis notebook analyzes the synthetic PROACT clinical dataset and the accompanying genetic variants (VCF).")
    add_code("""!pip install pandas numpy matplotlib seaborn statsmodels lifelines scikit-allel -q""")

    # Imports
    add_markdown("## 1. Import Libraries")
    add_code("""import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols
from lifelines import KaplanMeierFitter

# Set plot style
sns.set_theme(style="whitegrid")""")

    # Data Loading
    add_markdown("## 2. Load Clinical and Genetic Data")
    add_code("""# Load Clinical Data
clinical_df = pd.read_csv('synthetic_proact_dataset.csv')
print("Clinical Data Shape:", clinical_df.shape)
display(clinical_df.head())""")

    add_code("""# Load and parse VCF Data manually for simplicity
vcf_file = 'synthetic_cohort.vcf'
vcf_data = []

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            columns = line.strip().split('\\t')
            samples = columns[9:]
        else:
            parts = line.strip().split('\\t')
            variant_id = parts[2]
            genotypes = parts[9:]
            # Convert '0/0', '0/1', '1/1' to 0, 1, 2
            gt_numeric = [gt.count('1') for gt in genotypes]
            vcf_data.append({
                'Variant_ID': variant_id,
                'Genotypes': gt_numeric
            })

# Create Genotype DataFrame
genetic_df = pd.DataFrame({d['Variant_ID']: d['Genotypes'] for d in vcf_data})
genetic_df['Patient_ID'] = samples

# Merge Datasets
merged_df = pd.merge(clinical_df, genetic_df, on='Patient_ID')
print("Merged Data Shape:", merged_df.shape)
display(merged_df.head())""")

    # Clinical Analysis
    add_markdown("## 3. Clinical Data Analysis")
    add_code("""# Descriptive Statistics
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sns.histplot(data=merged_df, x='Age_at_Onset', kde=True, ax=axes[0])
axes[0].set_title('Age at Onset Distribution')

sns.countplot(data=merged_df, x='Site_of_Onset', ax=axes[1])
axes[1].set_title('Site of Onset')

sns.boxplot(data=merged_df, x='Site_of_Onset', y='ALSFRS_R_Slope', ax=axes[2])
axes[2].set_title('Progression Slope by Onset Site')

plt.tight_layout()
plt.show()""")

    add_code("""# Basic Survival Analysis (Bulbar vs Limb)
kmf = KaplanMeierFitter()
fig, ax = plt.subplots(figsize=(8, 5))

for site in merged_df['Site_of_Onset'].unique():
    mask = merged_df['Site_of_Onset'] == site
    kmf.fit(merged_df['Time_to_Event_Days'][mask], 
            event_observed=merged_df['Event_Occurred'][mask], 
            label=site)
    kmf.plot_survival_function(ax=ax)

plt.title('Kaplan-Meier Survival Curve by Site of Onset')
plt.ylabel('Survival Probability')
plt.xlabel('Days')
plt.show()""")

    # Genetic Subtyping
    add_markdown("## 4. Genetic Subtype Analysis\nWe define subtypes based on the Mendelian variants: `C9orf72` and `SOD1`.")
    add_code("""# Create Subtype Column
def classify_subtype(row):
    if row['C9orf72_exp'] > 0:
        return 'C9orf72'
    elif row['SOD1_A4V'] > 0 or row['SOD1_D90A'] > 0:
        return 'SOD1'
    else:
        return 'Sporadic'

merged_df['Genetic_Subtype'] = merged_df.apply(classify_subtype, axis=1)

print("Subtype Counts:")
print(merged_df['Genetic_Subtype'].value_counts())""")

    add_code("""# Phenotype verification: Site of Onset across Genetic Subtypes
# Based on our constraints: C9orf72 -> Bulbar, SOD1 -> Limb
plt.figure(figsize=(8, 5))
sns.countplot(data=merged_df, x='Genetic_Subtype', hue='Site_of_Onset')
plt.title('Site of Onset by Genetic Subtype')
plt.show()""")

    add_code("""# Survival Analysis by Genetic Subtype
fig, ax = plt.subplots(figsize=(8, 5))

for subtype in merged_df['Genetic_Subtype'].unique():
    mask = merged_df['Genetic_Subtype'] == subtype
    kmf.fit(merged_df['Time_to_Event_Days'][mask], 
            event_observed=merged_df['Event_Occurred'][mask], 
            label=subtype)
    kmf.plot_survival_function(ax=ax)

plt.title('Kaplan-Meier Survival Curve by Genetic Subtype')
plt.ylabel('Survival Probability')
plt.xlabel('Days')
plt.show()""")

    # GWAS Analysis
    add_markdown("## 5. GWAS Analysis on Sporadic Cohort\nTesting the polygenic risk alleles against progression slope (SARM1 constraint) and looking at overall frequencies.")
    add_code("""# Filter for Sporadic Cohort
sporadic_df = merged_df[merged_df['Genetic_Subtype'] == 'Sporadic'].copy()
print(f"Sporadic cohort size: {len(sporadic_df)}")

# Frequencies of Polygenic Risk Alleles
risk_alleles = ['rs12486783', 'rs11674437', 'rs11736730', 'rs35714649']

freqs = []
for snp in risk_alleles:
    # Allele frequency = (count of 1s + 2 * count of 2s) / (2 * N)
    af = (sporadic_df[snp].sum()) / (2 * len(sporadic_df))
    freqs.append({'SNP': snp, 'Allele Frequency': af})

freq_df = pd.DataFrame(freqs)
display(freq_df)""")

    add_code("""# Linear Regression: Testing SARM1 (rs35714649) association with Rapid Progression (ALSFRS_R_Slope)
# Since we injected this to be highly correlated with rapid progression (more negative slope)
# We expect a significant negative coefficient

# Add intercept for statsmodels
sporadic_df['Intercept'] = 1

X = sporadic_df[['Intercept', 'rs35714649', 'Age_at_Onset']]
y = sporadic_df['ALSFRS_R_Slope']

model = sm.OLS(y, X)
results = model.fit()

print(results.summary())

# Visualizing the association
plt.figure(figsize=(6, 5))
sns.boxplot(data=sporadic_df, x='rs35714649', y='ALSFRS_R_Slope')
plt.title('ALSFRS-R Slope by SARM1 (rs35714649) Genotype')
plt.xlabel('Genotype (0: Ref/Ref, 1: Ref/Alt)')
plt.ylabel('ALSFRS-R Slope (points/month)')
plt.show()""")

    # Write notebook
    with open("als_genotype_phenotype_analysis.ipynb", "w") as f:
        json.dump(notebook, f, indent=2)
    print("Notebook 'als_genotype_phenotype_analysis.ipynb' created successfully.")

if __name__ == "__main__":
    create_notebook()
