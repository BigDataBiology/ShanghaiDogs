import pandas as pd
from scipy.stats import shapiro,wilcoxon,friedmanchisquare
from itertools import combinations
from statsmodels.stats.multitest import multipletests


# Import files
polish_eval = pd.read_csv("intermediate-outputs/polishing_evaluation/results/Checkm2_final_polish_eval.csv", \
                  delimiter=',',header=0,index_col=0)
polish_eval_low = polish_eval[polish_eval['Mean_cov LR'] < 40]
polish_eval_high = polish_eval[polish_eval['Mean_cov LR'] >= 40]

# Assess normality - reshape to long format + Shapiro test
var = 'Completeness' #'Contamination' 'Completeness'
polish_df = polish_eval_low #polish_eval_high polish_eval_low polish_eval

df = polish_df[['Flye_completeness','Medaka_completeness','Polypolish_completeness','Polca_completeness']]
#df = polish_df[['Flye-contam','Medaka-contam','Poly-contam','Polca-contam']]
df = df.reset_index()
df = df.melt(id_vars="index", var_name="Condition", value_name=var)

conditions = df["Condition"].unique()
normality_results = {}
for condition in conditions:
    stat, p = shapiro(df[df["Condition"] == condition][var])
    normality_results[condition] = p
    print(f"Shapiro-Wilk test for {condition}: p-value = {p:.4f}")

### Non-normally distributed completeness, and contamination values

### Compute Friedmann + Wilcoxon test if significant
polish_df = polish_eval_low #polish_eval_high polish_eval_low polish_eval

# Perform Friedman test
stat, p = friedmanchisquare(polish_df['Flye_completeness'],
                                  polish_df['Medaka_completeness'],
                                  polish_df['Polypolish_completeness'],
                                  polish_df['Polca_completeness'])

print(f"Friedman test statistic = {stat:.4f}, p-value = {p:.4f}")

# Generate all consecutive pairs
#pairs = list(combinations(['Polca_completeness', 'Polypolish_completeness',
#                           'Medaka_completeness', 'Flye_completeness'], 2)) # Completeness
pairs = list(combinations(['Polca-contam', 'Poly-contam',
                           'Medaka-contam', 'Flye-contam'], 2)) # Contamination
print(pairs)

# Perform Wilcoxon for each pair
p_values = []
results = []

for step1, step2 in pairs:
    stat, p_value = wilcoxon(polish_eval_low[step1],polish_eval_low[step2],alternative='less') #'greater' for completeness
    results.append({'Comparison': f"{step1} vs {step2}", 'Statistic': stat, 'p-value': p_value, 'Group': 'low'})
    p_values.append(p_value)
    stat, p_value = wilcoxon(polish_eval_high[step1], polish_eval_high[step2],alternative='less') #'greater' for completeness
    results.append({'Comparison': f"{step1} vs {step2}", 'Statistic': stat, 'p-value': p_value, 'Group': 'high'})
    p_values.append(p_value)

# Correct for multiple testing (e.g., Benjamini-Hochberg)
_, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# Add corrected p-values to the results DataFrame
for idx, result in enumerate(results):
    result['Corrected p-value'] = p_values_corrected[idx]

# Create DataFrame
comparison_results = pd.DataFrame(results)
comparison_results['Corrected p-value'] = comparison_results['Corrected p-value'].apply(lambda x: f"{x:.2e}")
