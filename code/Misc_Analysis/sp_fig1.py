import pandas as pd
import matplotlib.pyplot as plt

stats_data = pd.read_csv('/mnt/d/hanjia/2.10/new_fig5/time_diff_stats_with_filtered.tsv', sep='\t', index_col=0)

fig, ax = plt.subplots(figsize=(6, 4))

metrics = ['Mean', 'Median', 'Q1', 'Q3']
values = [stats_data.loc['Mean', 'time_diff_days'], 
          stats_data.loc['Median', 'time_diff_days'],
          stats_data.loc['Q1', 'time_diff_days'],
          stats_data.loc['Q3', 'time_diff_days']]

colors = ['#2ecc71', '#3498db', '#9b59b6', '#e74c3c']

values = [-x for x in values]
bars = ax.barh(metrics, values, color=colors, alpha=0.7)

ax.set_title('Time Difference Statistics', fontsize=12)
ax.set_xlabel('Days', fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.grid(axis='x', linestyle='--', alpha=0.3)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

for bar in bars:
    width = bar.get_width()
    ax.text(0, bar.get_y() + bar.get_height()/2.,
            f'{-width:.1f}', ha='left', va='center', fontsize=10)

plt.tight_layout()

plt.savefig('/mnt/d/hanjia/2.10/new_fig5/sp_fig1.pdf')
plt.savefig('/mnt/d/hanjia/2.10/new_fig5/sp_fig1.svg')
plt.show()