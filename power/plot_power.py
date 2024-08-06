#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#%%
data = pd.read_csv("vary_genes_incld_self_pval.csv", delimiter=' ')
row_counts = data['h2'].value_counts().reset_index()
row_counts.columns = ['h2', 'row_count']
#%%
data
# %%
threshold1 = 0.05/9515
result1 = data.groupby('h2')['pval'].apply(lambda x: (x < threshold1).mean()).reset_index()
result1.columns = ['h2', 'pval_ratio']

error_vals1 = []
for row, pair in result1.iterrows():
    h2_counts = row_counts.loc[row_counts['h2'] == pair['h2'], 'row_count'].values[0]
    error = 2 * np.sqrt(pair['pval_ratio'] * (1 - pair['pval_ratio']) / h2_counts)
    error_vals1.append(error)
print(error_vals1)
#%%
threshold2 = 0.05
result2 = data.groupby('h2')['pval'].apply(lambda x: (x < threshold2).mean()).reset_index()
result2.columns = ['h2', 'pval_ratio']

error_vals2 = []
for row, pair in result2.iterrows():
    h2_counts = row_counts.loc[row_counts['h2'] == pair['h2'], 'row_count'].values[0]
    error = 2 * np.sqrt(pair['pval_ratio'] * (1 - pair['pval_ratio']) / h2_counts)
    error_vals2.append(error)
print(error_vals2)
# %%
f, ax = plt.subplots(figsize=(8, 10))
ax.tick_params(axis='x', which='major')
plt.errorbar([0, 1, 2, 3, 4], result1['pval_ratio'], yerr=error_vals1, linestyle='-', marker='o', markersize=10, linewidth=2, capsize=10, label=r'${\alpha = \frac{0.05}{9515}}$', color='red')
plt.errorbar([0, 1, 2, 3, 4], result2['pval_ratio'], yerr=error_vals2, linestyle='-', marker='^', markersize=10, linewidth=2, capsize=10, label=r'${\alpha = 0.05}$', color='blue')
plt.xticks(ticks=[0, 1, 2, 3, 4], labels=[0.0001, 0.0005, 0.001, 0.005, 0.01])
plt.ylim([0, 1.05])

plt.ylabel("Power", size = 30)
plt.xlabel("$h^{2}_{quad}$", size = 30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=30,markerscale=1.,loc='lower right')
# plt.text(weight=1)
#plt.subplots_adjust(bottom=0.18)
#plt.subplots_adjust(left=0.15)
plt.savefig("self_interaction_incld.png", dpi=200, bbox_inches='tight')
# %%
