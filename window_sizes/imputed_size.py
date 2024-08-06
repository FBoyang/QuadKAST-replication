#%%
import pandas as pd
import numpy as np
#%%
array_data = pd.read_csv("/u/home/a/aanand2/window_size_testing/genes_info_array")
imputed_data = pd.read_csv("/u/home/a/aanand2/window_size_testing/qced.bim", delimiter="\t", header=None)
# %%
chrs = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 3, 4, 5, 6, 7, 8, 9]
chr_ind = 0
chr = chrs[chr_ind]
prev_start = 0

inds = array_data.loc[:, ["pStart","pEnd"]]
imputed_window_sizes = []
genes = []
g = 0

for i,j in inds.values:
    print(g)
    gene = array_data.loc[g, "gene"]
    genes.append(gene)
    start = i
    end = j

    if start < prev_start:
        chr_ind += 1
        chr = chrs[chr_ind]

    wsize = sum(imputed_data[imputed_data[0]==int(chr)][3].between(int(start), int(end)))
    imputed_window_sizes.append(wsize)

    prev_start = start

    g += 1
# %%
df = pd.DataFrame({"gene":genes,"#SNPs":imputed_window_sizes})
df.to_csv("/u/home/a/aanand2/window_size_testing/imputed_sizes", index=False)
#np.savetxt("/u/home/a/aanand2/window_size_testing/imputed_sizes.txt", imputed_window_sizes)
# %%
