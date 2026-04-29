import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
import base64
from io import BytesIO
import os

# Set style
sns.set_theme(style="whitegrid")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Inter', 'Arial']

def fig_to_base64(fig):
    img = BytesIO()
    fig.savefig(img, format='png', bbox_inches='tight', dpi=150)
    img.seek(0)
    return base64.b64encode(img.getvalue()).decode('utf-8')

def generate_report():
    # Load data
    csv_path = 'synthetic_proact_dataset.csv'
    df = pd.read_csv(csv_path)
    
    # 1. Age of Onset Distribution
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    sns.histplot(df['Age_at_Onset'], kde=True, color='#6366f1', ax=ax1)
    ax1.set_title('Distribution of Age of Onset', fontsize=16, fontweight='bold', pad=20)
    ax1.set_xlabel('Age (years)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    age_img = fig_to_base64(fig1)
    plt.close(fig1)

    # 2. Site of Onset Distribution
    fig2, ax2 = plt.subplots(figsize=(8, 8))
    site_counts = df['Site_of_Onset'].value_counts()
    colors = ['#818cf8', '#fb7185', '#34d399', '#fbbf24']
    ax2.pie(site_counts, labels=site_counts.index, autopct='%1.1f%%', 
            colors=colors, startangle=140, textprops={'fontsize': 12})
    ax2.set_title('Distribution of Site of Onset', fontsize=16, fontweight='bold', pad=20)
    site_img = fig_to_base64(fig2)
    plt.close(fig2)

    # 3. Progression by Site of Onset
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    sns.boxplot(x='Site_of_Onset', y='ALSFRS_R_Slope', data=df, palette=colors, ax=ax3)
    ax3.set_title('Disease Progression (ALSFRS-R Slope) by Site of Onset', fontsize=16, fontweight='bold', pad=20)
    ax3.set_xlabel('Site of Onset', fontsize=12)
    ax3.set_ylabel('Slope (points/month)', fontsize=12)
    prog_img = fig_to_base64(fig3)
    plt.close(fig3)

    # 4. Survival by Site of Onset (Kaplan-Meier)
    fig4, ax4 = plt.subplots(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    for i, site in enumerate(df['Site_of_Onset'].unique()):
        mask = df['Site_of_Onset'] == site
        kmf.fit(df.loc[mask, 'Time_to_Event_Days'], 
                event_observed=df.loc[mask, 'Event_Occurred'], 
                label=site)
        kmf.plot_survival_function(ax=ax4, color=colors[i % len(colors)])
        
    ax4.set_title('Survival by Site of Onset (Kaplan-Meier Curves)', fontsize=16, fontweight='bold', pad=20)
    ax4.set_xlabel('Time (Days)', fontsize=12)
    ax4.set_ylabel('Survival Probability', fontsize=12)
    ax4.legend(title='Site of Onset')
    survival_img = fig_to_base64(fig4)
    plt.close(fig4)

    # Statistics Summary
    stats = {
        'total_patients': len(df),
        'mean_age': df['Age_at_Onset'].mean(),
        'median_survival': df['Time_to_Event_Days'].median(),
        'site_breakdown': df['Site_of_Onset'].value_counts().to_dict()
    }

    # HTML Template
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ALS Clinical Data Analysis Report</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --primary: #6366f1;
            --text-main: #1e293b;
            --text-muted: #64748b;
            --border: #e2e8f0;
        }}
        
        body {{
            font-family: 'Inter', sans-serif;
            background-color: var(--bg-color);
            color: var(--text-main);
            margin: 0;
            padding: 40px 20px;
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1000px;
            margin: 0 auto;
        }}
        
        header {{
            text-align: center;
            margin-bottom: 60px;
        }}
        
        h1 {{
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 10px;
            color: var(--primary);
        }}
        
        .subtitle {{
            font-size: 1.1rem;
            color: var(--text-muted);
        }}
        
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(450px, 1fr));
            gap: 30px;
            margin-bottom: 40px;
        }}
        
        .card {{
            background: var(--card-bg);
            border-radius: 16px;
            padding: 30px;
            box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1), 0 2px 4px -2px rgb(0 0 0 / 0.1);
            border: 1px solid var(--border);
            transition: transform 0.2s;
        }}
        
        .card:hover {{
            transform: translateY(-4px);
        }}
        
        .card h2 {{
            font-size: 1.25rem;
            font-weight: 600;
            margin-bottom: 20px;
            border-left: 4px solid var(--primary);
            padding-left: 15px;
        }}
        
        img {{
            width: 100%;
            height: auto;
            border-radius: 8px;
        }}
        
        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 20px;
            margin-bottom: 40px;
        }}
        
        .stat-box {{
            background: var(--primary);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
        }}
        
        .stat-value {{
            font-size: 2rem;
            font-weight: 700;
            display: block;
        }}
        
        .stat-label {{
            font-size: 0.875rem;
            opacity: 0.9;
        }}

        footer {{
            text-align: center;
            margin-top: 60px;
            color: var(--text-muted);
            font-size: 0.875rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>ALS Clinical Data Analysis</h1>
            <p class="subtitle">Comprehensive phenotypic report for cohort PRO-ACT-SYN</p>
        </header>
        
        <div class="summary-stats">
            <div class="stat-box">
                <span class="stat-value">{stats['total_patients']}</span>
                <span class="stat-label">Total Patients</span>
            </div>
            <div class="stat-box" style="background: #ec4899;">
                <span class="stat-value">{stats['mean_age']:.1f}</span>
                <span class="stat-label">Mean Age at Onset</span>
            </div>
            <div class="stat-box" style="background: #10b981;">
                <span class="stat-value">{stats['median_survival']:.0f}</span>
                <span class="stat-label">Median Survival (Days)</span>
            </div>
        </div>
        
        <div class="grid">
            <div class="card">
                <h2>Age at Onset</h2>
                <img src="data:image/png;base64,{age_img}" alt="Age Distribution">
                <p>The distribution shows a peak onset around {stats['mean_age']:.1f} years, consistent with typical ALS demographics.</p>
            </div>
            
            <div class="card">
                <h2>Site of Onset</h2>
                <img src="data:image/png;base64,{site_img}" alt="Site of Onset Distribution">
                <p>Breakdown of onset location across the cohort, showing the relative frequency of bulbar vs limb onset.</p>
            </div>
            
            <div class="card">
                <h2>Disease Progression</h2>
                <img src="data:image/png;base64,{prog_img}" alt="Progression Slope">
                <p>Stratified ALSFRS-R slope analysis showing functional decline across different onset sites.</p>
            </div>
            
            <div class="card">
                <h2>Survival Analysis</h2>
                <img src="data:image/png;base64,{survival_img}" alt="Kaplan-Meier Survival">
                <p>Kaplan-Meier curves illustrating survival probability differences between onset sites.</p>
            </div>
        </div>
        
        <footer>
            <p>Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')} | Data Source: synthetic_proact_dataset.csv</p>
        </footer>
    </div>
</body>
</html>
    """

    with open('als_analysis_report.html', 'w') as f:
        f.write(html_content)
    
    print("Report generated successfully: als_analysis_report.html")

if __name__ == "__main__":
    generate_report()
