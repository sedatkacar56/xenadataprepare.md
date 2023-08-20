```import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statannot import add_stat_annotation```


# Read TSV data with headers
data = pd.read_csv("C:/Users/sedat/Downloads/ccl_GTEX/ccl2.tsv", sep='\t')

# Define the sample types to include
sample_types_to_include = ["Primary Tumor", "Solid Tissue Normal", "Normal Tissue"]

# Filter data and select relevant columns
filtered_data = data[data['_sample_type'].isin(sample_types_to_include)][['_sample_type', data.columns[2]]]

# Rename columns for clarity
filtered_data.columns = ['Sample_Type', 'Gene_Expression']  # Change 'Gene_Expression' to the correct column name

# Define the desired order of categories
desired_order = ["Normal Tissue", "Solid Tissue Normal", "Primary Tumor"]

# Perform ANOVA
groups = filtered_data['Sample_Type'].unique()
anova_results = {}
for group in groups:
    group_data = filtered_data[filtered_data['Sample_Type'] == group]['Gene_Expression']
    anova_results[group] = group_data

f_statistic, p_value = f_oneway(*anova_results.values())
print("F-statistic:", f_statistic)
print("P-value:", p_value)

# Perform Tukey's HSD post hoc test
posthoc = pairwise_tukeyhsd(filtered_data['Gene_Expression'], filtered_data['Sample_Type'])
print(posthoc)

# Create box plot and add asterisks for significant differences
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
ax = sns.boxplot(x='Sample_Type', y='Gene_Expression', data=filtered_data, order=desired_order, showfliers=False)
#plt.title("Box Plot of Gene Expression by Sample Type")
plt.xlabel("Sample Type")
plt.ylabel("Gene Expression")



plt.show()
