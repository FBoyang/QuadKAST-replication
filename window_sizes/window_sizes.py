#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%%
array_data = pd.read_csv("/u/home/a/aanand2/window_size_testing/genes_info_array")
# %%
array_data_filt_genes = array_data[array_data["#SNPs"] <= 50].loc[:, "gene"]
len(array_data_filt_genes.index)
#%%
window_sizes = array_data.loc[:, "#SNPs"]
percentile = np.percentile(window_sizes, 99)
print("99th percentile: " + str(percentile))
print("Median: " + str(np.median(window_sizes)))
print(len(window_sizes))
window_sizes = window_sizes[window_sizes < percentile]
print(len(window_sizes))
# %%
plt.figure(figsize=(8, 6), dpi=80)
plt.hist(window_sizes, bins=40, weights=np.ones(len(window_sizes)) / len(window_sizes)*100)
plt.axvline(50, color='r', linestyle='dashed', linewidth=1)
plt.xlabel("Window Size (SNPs)", size=20)
plt.ylabel("% Genes", size=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(prop={'size': 16},markerscale=2.,loc='upper right')
plt.savefig('/u/home/a/aanand2/window_size_testing/' + "array_size" + '.png', bbox_inches='tight', dpi=200)
#%%
imputed_window_sizes_list = np.loadtxt("/u/home/a/aanand2/window_size_testing/imputed_sizes.txt")
imputed_window_sizes = pd.DataFrame(imputed_window_sizes_list, columns=["#SNPs"])
print(imputed_window_sizes)
#%%
imputed_window_sizes_filt_genes = [imputed_window_sizes_list[i] for i in array_data_filt_genes.index.values.tolist()]
len(imputed_window_sizes_filt_genes)
#%%
imputed_window_sizes_filt_genes.index(max(imputed_window_sizes_filt_genes))
#%%
percentile = np.percentile(imputed_window_sizes, 99)
print("percentile: " + str(percentile))
print("Median: " + str(np.median(imputed_window_sizes)))

imputed_window_sizes = imputed_window_sizes.loc[imputed_window_sizes["#SNPs"] <= percentile]
#print(imputed_window_sizes)
# %%
plt.figure(figsize=(8, 6), dpi=80)
plt.hist(imputed_window_sizes, bins=40, weights=np.ones(len(imputed_window_sizes)) / len(imputed_window_sizes)*100)
#plt.axvline(1542, color='r', linestyle='dashed', linewidth=1)
#plt.hist(imputed_window_sizes, bins=40)
plt.xlabel("Window Size (SNPs)", size=20)
plt.ylabel("% Genes", size=20)
#plt.title("Imputed Data", size=25)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(prop={'size': 16},markerscale=2.,loc='upper right')
# %%
exome_window_sizes = np.loadtxt("/u/home/a/aanand2/window_size_testing/exome_sizes.txt")
exome_window_sizes = pd.DataFrame(exome_window_sizes, columns=["#SNPs"])
exome_window_sizes = exome_window_sizes.loc[exome_window_sizes["#SNPs"] != 0]
exome_window_sizes = exome_window_sizes.loc[exome_window_sizes["#SNPs"] < 6000]
print(len(exome_window_sizes))
# %%
plt.figure(figsize=(8, 6), dpi=80)
plt.hist(exome_window_sizes, bins=40, weights=np.ones(len(exome_window_sizes)) / len(exome_window_sizes)*100)
plt.xlabel("Window Size (SNPs)", size=20)
plt.ylabel("% Genes", size=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(prop={'size': 16},markerscale=2.,loc='upper right')
# %%
