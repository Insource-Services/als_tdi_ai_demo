# Genetic Subtype Analysis Plan

This plan details the addition of the genetic subtype analysis to the `als_genotype_phenotype_analysis.ipynb` notebook. The goal is to evaluate the differences in site of onset and survival probability among patients with SOD1, C9orf72, and sporadic ALS.

## Open Questions

- If a patient has both C9orf72 and SOD1 variants (which is biologically rare but possible in synthetic data), they will be labeled according to the check order. I propose prioritizing C9orf72, followed by SOD1. Does this work for you?
- The current dataframe has `SOD1_A4V` and `SOD1_D90A` columns. I will define a patient as `SOD1` if either of these variants are > 0. Does this cover the SOD1 genotype definition correctly for this synthetic dataset?

## Proposed Changes

### Notebook Modifications

#### [MODIFY] [als_genotype_phenotype_analysis.ipynb](file:///Users/bcole/bin/als_tdi_ai_demo/als_genotype_phenotype_analysis.ipynb)
- **Cell Modification 1 (Data Prep & Count Plot)**: Replace the first empty code cell under the `## 4. Genetic Subtype Analysis` header with Python code to:
  - Create a new `Genotype` column based on `merged_df` conditions: 
    - `C9orf72` if `C9orf72_exp > 0`
    - `SOD1` if `SOD1_A4V > 0` or `SOD1_D90A > 0`
    - `Sporadic` otherwise.
  - Create a bar plot using `sns.countplot(data=merged_df, x='Site_of_Onset', hue='Genotype')`.
- **Cell Modification 2 (Survival Analysis)**: Replace the second empty code cell with Python code to:
  - Initialize `KaplanMeierFitter`.
  - Loop over the unique values of `Genotype` (`C9orf72`, `SOD1`, `Sporadic`).
  - Fit and plot the survival function based on `Time_to_Event_Days` and `Event_Occurred` for each genotype to display their respective survival curves on a single plot.

## Verification Plan

### Manual Verification
- After inserting the analysis code, I can execute the notebook using `jupyter nbconvert` or just let you run the cells in your environment to view the generated plots and verify the data distribution.
