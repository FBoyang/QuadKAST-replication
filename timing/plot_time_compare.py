#%%
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
#%%
M = 100
N = [0, 5000, 10000, 20000, 50000, 100000, 200000, 291273]
Kernels = ["SKAT_QUAD", "DIRECT"]
Windows = [100]

sample_sizes = []
window_sizes = []
kernel = []
times = []
mems = []

for k in Kernels:
    for w in Windows:
        sample_sizes.append(0)
        window_sizes.append(w)
        kernel.append(k)
        times.append(0)
        mems.append(0)

basedir = f"/u/home/a/aanand2/QuadKAST-replication/timing/{M}snp/"
trials = 10
for k in Kernels:
    for t in range(1, trials+1):
        timedir = basedir + "time/" + k + "/t" + str(t) + "/"
        print(timedir)
        if (not os.path.isdir(timedir)):
            continue
        for filename in os.listdir(timedir):
            f = os.path.join(timedir, filename)
            timefile = open(f)
            time = float(timefile.readline().strip())/3600
            timefile.close()
            memfile = open(basedir + "mem/" + k + "/t" + str(t) + "/" + filename)
            content = memfile.read().splitlines()
            content = [c for c in content if c]
            startmem = content[content.index("=============================================================")+1].split()[1]
            endmem = content[len(content)-1].split()[1]
            mem = float(endmem) - float(startmem)
            memfile.close()

            sample_sizes.append(N[int(filename[0])])
            window_sizes.append(w)
            kernel.append(k)
            times.append(time)
            mems.append(mem)

d = {"Sample Size" : sample_sizes, "Window Size" : window_sizes, "Kernel" : kernel, "Time (hrs)" : times, "Memory (MiB)" : mems}
df = pd.DataFrame(d)

time_df = df.groupby(["Sample Size", "Kernel"])["Time (hrs)"].mean().reset_index(name='Time (hrs)')
time_df_std = df.groupby(["Sample Size", "Kernel"])["Time (hrs)"].std().reset_index(name='std')
time_df_std["std"] = time_df_std["std"]/np.sqrt(10)
time_df["std"] = time_df_std["std"]
time_df.fillna(0, inplace=True)

mem_df = df.groupby(["Sample Size", "Kernel"])["Memory (MiB)"].mean().reset_index(name='Memory (MiB)')
mem_df_std = df.groupby(["Sample Size", "Kernel"])["Memory (MiB)"].std().reset_index(name='std')
mem_df_std["std"] = mem_df_std["std"]/np.sqrt(10)
mem_df["std"] = mem_df_std["std"]
mem_df.fillna(0, inplace=True)
# print(time_df)
# print(mem_df)
#%%
f, ax = plt.subplots(figsize=(8, 10), facecolor="w", edgecolor="xkcd:grey")
plt.plot(time_df[time_df["Kernel"] == "SKAT_QUAD"]["Sample Size"]/1000, time_df[time_df["Kernel"] == "SKAT_QUAD"]["Time (hrs)"], label="SKAT_QUAD", marker='o', markersize=10, linewidth=4, color='blue')
plt.plot(time_df[time_df["Kernel"] == "DIRECT"]["Sample Size"]/1000, time_df[time_df["Kernel"] == "DIRECT"]["Time (hrs)"], label="QuadKAST", marker='^', markersize=10, linewidth=4, color='red')
#plt.xlabel("Sample Size (k)", size=14)
#plt.ylabel("Time (hrs)", size=14)
#plt.title("Time (M = 100)", size=25)
plt.xticks(fontsize=30, weight=750)
plt.yticks(fontsize=30, weight=750)
plt.legend(prop={'size': 30, 'weight' : 750},markerscale=2.,loc='upper right')
plt.subplots_adjust(bottom=0.17)
#plt.subplots_adjust(left=0.15)
# plt.savefig("plots/time_compr_m100.png", dpi=200)
#%%
f, ax = plt.subplots(figsize=(6, 3.5), facecolor="w", edgecolor="k")
plt.plot(mem_df[mem_df["Kernel"] == "SKAT_QUAD"]["Sample Size"]/1000, mem_df[mem_df["Kernel"] == "SKAT_QUAD"]["Memory (MiB)"]/953.7, label="SKAT_QUAD", marker='o')
plt.plot(mem_df[mem_df["Kernel"] == "DIRECT"]["Sample Size"]/1000, mem_df[mem_df["Kernel"] == "DIRECT"]["Memory (MiB)"]/953.7, label="QuadKAST", marker='^')
plt.xlabel("Sample Size (k)", size=14)
plt.ylabel("Memory (GB)", size=14)
#plt.title("Memory (M = 100)", size=25)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 13},markerscale=2.,loc='upper right')
plt.subplots_adjust(bottom=0.17)
#plt.subplots_adjust(left=0.15)
# plt.savefig("plots/mem_compr_m100.png", dpi=200)
#%%
new_df = df[df["Kernel"]=="DIRECT"][1:]

new_df['Time (hrs)'] = new_df['Time (hrs)'] * 60
# Group by 'Sample Size' and collect 'Time (hrs)' into lists
grouped = new_df.groupby('Sample Size')['Time (hrs)'].apply(list).reset_index()

# Expand lists into separate columns (Trial 1 to Trial 10)
expanded = pd.DataFrame(grouped['Time (hrs)'].tolist(), index=grouped['Sample Size'])
expanded.columns = [f'Trial {i+1}' for i in range(expanded.shape[1])]
#%%
expanded.to_csv("/u/home/a/aanand2/QuadKAST-replication/timing/runtimes_sample.csv")
# %%
## SNP TEST ##
# %%
import os
import pandas as pd
#%%
# Define the directory containing the data
data_directory = '/u/home/a/aanand2/QuadKAST-replication/timing/50000n'

# Initialize an empty DataFrame to store the results
results_df = pd.DataFrame(columns=[10, 50, 100, 200, 300, 500])
#%%
# Loop through the trial folders (t1 to t10)
for trial_folder in range(1, 11):
    trial_folder_path = os.path.join(data_directory, f't{trial_folder}')
    
    # Initialize a list to store the results for the current trial
    trial_results = []
    
    # Loop through the setting files (1.txt to 6.txt)
    for setting_file in range(1, 7):
        setting_file_path = os.path.join(trial_folder_path, f'{setting_file}.txt')
        
        # Read the number from the file and append it to the trial results list
        with open(setting_file_path, 'r') as file:
            number = float(file.read().strip())
            trial_results.append(number/60)
    
    # Append the trial results to the results DataFrame with the trial number as the index
    results_df.loc[f'Trial {trial_folder}'] = trial_results

# Transpose the DataFrame to have the settings as columns and trials as rows
results_df = results_df.transpose()
#%%
# Rename the index and reset it to start from 10
results_df.index.name = 'Window Size'
results_df.reset_index(inplace=True)
results_df['Window Size'] = [10, 50, 100, 200, 300, 500]
# %%
results_df.to_csv("runtimes_snp.csv", index=False)
# %%