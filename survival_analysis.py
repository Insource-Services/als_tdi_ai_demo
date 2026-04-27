import pandas as pd
import numpy as np

def generate_proact_synthetic_data(n_patients=1000, random_seed=42):
    """
    Generates a synthetic dataset of ALS patients based on characteristics 
    from the PRO-ACT (Pooled Resource Open-Access ALS Clinical Trials) database.
    """
    np.random.seed(random_seed)
    
    # 1. Demographics & Baseline Characteristics
    # Age at onset: mean ~56, std ~12 (normal distribution, bounded)
    age_onset = np.random.normal(loc=56, scale=12, size=n_patients)
    age_onset = np.clip(age_onset, 20, 90) # Clip to realistic ranges
    
    # Sex: ~40% Female, 60% Male
    sex = np.random.choice(['Male', 'Female'], size=n_patients, p=[0.60, 0.40])
    
    # Site of Onset: ~76% Limb, 24% Bulbar
    onset_site = np.random.choice(['Limb', 'Bulbar'], size=n_patients, p=[0.76, 0.24])
    
    # Riluzole Use: ~78% Yes, 22% No
    riluzole_use = np.random.choice(['Yes', 'No'], size=n_patients, p=[0.78, 0.22])
    
    # 2. Clinical Measurements
    # ALSFRS-R Baseline: max 48, typical mean ~38, std ~6
    # We'll use a beta distribution scaled to 0-48 to get a left-skewed distribution
    alsfrs_baseline = 48 - np.random.gamma(shape=2.5, scale=4, size=n_patients)
    alsfrs_baseline = np.clip(np.round(alsfrs_baseline), 10, 48)
    
    # Baseline FVC (Forced Vital Capacity % predicted): mean ~85, std ~20
    fvc_baseline = np.random.normal(loc=85, scale=20, size=n_patients)
    fvc_baseline = np.clip(np.round(fvc_baseline), 30, 130)
    
    # 3. Disease Progression & Associations
    # ALSFRS-R Slope (points/month): typical mean ~ -1.02, std ~0.8
    # Slope is associated with onset site (bulbar faster), age (older faster)
    
    base_slope = np.random.normal(loc=-0.8, scale=0.5, size=n_patients)
    
    # Add penalties (faster decline) for bulbar onset and older age
    onset_penalty = np.where(onset_site == 'Bulbar', -0.3, 0)
    age_penalty = (age_onset - 56) * -0.015  # older = more negative slope
    
    alsfrs_slope = base_slope + onset_penalty + age_penalty
    # Clip to realistic bounds
    alsfrs_slope = np.clip(alsfrs_slope, -4.0, 0.5) 
    
    # 4. Survival Outcomes
    # We'll use a Cox proportional hazards-like model to generate survival times
    # Median survival from trial entry is ~480 days, but let's model from symptom onset or trial entry
    # Let's generate survival time from trial baseline (days)
    
    # Calculate a risk score (log hazard) for each patient
    # Higher score = higher risk of death/event = shorter survival
    risk_score = np.zeros(n_patients)
    
    # Age effect (older = higher risk)
    risk_score += (age_onset - 56) * 0.03
    
    # Onset site effect (bulbar = higher risk)
    risk_score += np.where(onset_site == 'Bulbar', 0.5, 0)
    
    # Baseline ALSFRS effect (lower score = higher risk)
    risk_score += (38 - alsfrs_baseline) * 0.05
    
    # Baseline FVC effect (lower FVC = higher risk)
    risk_score += (85 - fvc_baseline) * 0.02
    
    # Slope effect (more negative slope = faster decline = higher risk)
    risk_score += (-1.0 - alsfrs_slope) * 0.8
    
    # Convert risk score to a hazard multiplier
    hazard_multiplier = np.exp(risk_score)
    
    # Generate baseline survival times from a Weibull distribution
    # shape parameter (k) > 1 indicates increasing hazard over time
    baseline_scale = 800  # determines overall median survival
    shape = 1.5
    
    # Inverse transform sampling for Weibull with proportional hazards
    u = np.random.uniform(0, 1, size=n_patients)
    survival_days = baseline_scale * (-np.log(u) / hazard_multiplier) ** (1 / shape)
    
    # 5. Censoring
    # Clinical trials typically have a fixed follow-up period (e.g., 12-18 months / 365-550 days)
    # plus some patients lost to follow-up
    max_follow_up = np.random.uniform(300, 700, size=n_patients)
    
    # An event occurs if survival time is less than follow-up time
    event_occurred = (survival_days <= max_follow_up).astype(int)
    
    # Observed time is the minimum of survival time and follow-up time
    time_to_event = np.minimum(survival_days, max_follow_up)
    
    # Compile into a DataFrame
    df = pd.DataFrame({
        'Patient_ID': [f'PROACT_SYN_{i:04d}' for i in range(1, n_patients + 1)],
        'Age_at_Onset': np.round(age_onset, 1),
        'Sex': sex,
        'Site_of_Onset': onset_site,
        'Riluzole_Use': riluzole_use,
        'ALSFRS_R_Baseline': alsfrs_baseline,
        'FVC_Baseline': fvc_baseline,
        'ALSFRS_R_Slope': np.round(alsfrs_slope, 3),
        'Time_to_Event_Days': np.round(time_to_event),
        'Event_Occurred': event_occurred
    })
    
    return df

if __name__ == "__main__":
    print("Generating synthetic PROACT dataset with 1,000 patients...")
    df_synthetic = generate_proact_synthetic_data(n_patients=1000)
    
    # Print some summary statistics to verify distributions
    print("\nDataset Shape:", df_synthetic.shape)
    print("\nSummary Statistics:")
    print(df_synthetic.describe())
    print("\nCategorical Distributions:")
    print(df_synthetic['Sex'].value_counts(normalize=True))
    print(df_synthetic['Site_of_Onset'].value_counts(normalize=True))
    
    # Save to CSV
    output_file = "synthetic_proact_dataset.csv"
    df_synthetic.to_csv(output_file, index=False)
    print(f"\nSynthetic dataset saved to {output_file}")
