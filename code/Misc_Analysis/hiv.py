import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20

df_2023_age = pd.read_csv('/mnt/d/hanjia/2.10/Baseline/2023_age.csv')
df_2024_age = pd.read_csv('/mnt/d/hanjia/2.10/Baseline/2024_age.csv')
df_2023_status = pd.read_csv('/mnt/d/hanjia/2.10/fig_4/2023.csv')
df_2024_status = pd.read_csv('/mnt/d/hanjia/2.10/fig_4/2024.csv')

df_combined_age = pd.concat([df_2023_age, df_2024_age])

fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))

ax1.boxplot([df_2023_age['age'], df_2024_age['age'], df_combined_age['age']], 
            labels=['2023', '2024', 'Combined'])
ax1.set_ylabel('Age')
# ax1.set_title('Age Distribution by Year')

plt.tight_layout()
plt.savefig('age_distribution.pdf', dpi=1000, bbox_inches='tight')
plt.savefig('age_distribution.svg', format='svg', bbox_inches='tight')
plt.savefig('age_distribution.png', dpi=1000, bbox_inches='tight')
plt.show()

fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))

total_2023 = len(df_2023_status)
total_2024 = len(df_2024_status)
hiv_2023 = len(df_2023_status[df_2023_status['HIV'] == 'Yes'])
hiv_2024 = len(df_2024_status[df_2024_status['HIV'] == 'Yes'])
msm_2024 = len(df_2024_status[df_2024_status['MSM'] == 'Yes'])

hiv_ratio_2023 = (hiv_2023 / total_2023) * 100
hiv_ratio_2024 = (hiv_2024 / total_2024) * 100
msm_ratio_2024 = (msm_2024 / total_2024) * 100

x = np.array([1, 2, 3])
width = 0.35

ax2.bar(x, [100, 100, 100], width, color='#DCE9F4', label='Total')
ax2.bar(x[0], hiv_ratio_2023, width, color='#43978F', label='2023 HIV')
ax2.bar(x[1], hiv_ratio_2024, width, color='#9EC4BE', label='2024 HIV')
ax2.bar(x[2], msm_ratio_2024, width, color='#ABD0F1', label='2024 MSM')

ax2.set_ylim(0, 100)
ax2.set_ylabel('Percentage (%)')
ax2.set_xticks(x)
ax2.set_xticklabels(['2023 HIV', '2024 HIV', '2024 MSM'])
# ax2.set_title('HIV and MSM Prevalence')
ax2.legend()

plt.tight_layout()
plt.savefig('hiv_msm_prevalence.pdf', dpi=1000, bbox_inches='tight')
plt.savefig('hiv_msm_prevalence.svg', format='svg', bbox_inches='tight')
plt.savefig('hiv_msm_prevalence.png', dpi=1000, bbox_inches='tight')
plt.show()

with open('separate_stats.txt', 'w') as f:
    f.write("=== Age Distribution Statistics ===\n")
    f.write(f"2023 Sample size: {len(df_2023_age)}\n")
    f.write(f"2024 Sample size: {len(df_2024_age)}\n")
    f.write(f"Combined sample size: {len(df_combined_age)}\n")
    f.write(f"2023 Average age: {df_2023_age['age'].mean():.1f}\n")
    f.write(f"2024 Average age: {df_2024_age['age'].mean():.1f}\n")
    f.write(f"Combined average age: {df_combined_age['age'].mean():.1f}\n\n")
    
    f.write("=== HIV/MSM Prevalence Statistics ===\n")
    f.write(f"2023 HIV positive prevalence: {hiv_ratio_2023:.1f}%\n")
    f.write(f"2024 HIV positive prevalence: {hiv_ratio_2024:.1f}%\n")
    f.write(f"2024 MSM prevalence: {msm_ratio_2024:.1f}%\n")

print("Plots have been saved as:")
print("- Left plot (Age distribution): age_distribution.pdf, age_distribution.svg")
print("- Right plot (HIV/MSM prevalence): hiv_msm_prevalence.pdf, hiv_msm_prevalence.svg")
print("- Statistics: separate_stats.txt")