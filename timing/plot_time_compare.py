#%%
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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


# def quad_line(x, a, b, c):
#     return a * x *x  + b*x + c

def quad_line(x, a, b, c, d, e):
    return a * x**4 + b*x**3 + c*x**2  + d*x + e

#%%
for k in Kernels:
    for w in Windows:
        sample_sizes.append(0)
        window_sizes.append(w)
        kernel.append(k)
        times.append(0)
        mems.append(0)

basedir = f"{M}snp/"
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
df["Time (mins)"]  = df["Time (hrs)"] * 60
time_df = df.groupby(["Sample Size", "Kernel"])["Time (mins)"].mean().reset_index(name='Time (mins)')
time_df_std = df.groupby(["Sample Size", "Kernel"])["Time (mins)"].std().reset_index(name='std')
time_df_std["std"] = time_df_std["std"]/np.sqrt(10)
time_df["std"] = time_df_std["std"]
time_df.fillna(0, inplace=True)

time_df_disp = time_df[time_df['Sample Size']>0].sort_values(by=['Kernel','Sample Size'])
time_df_disp.loc[time_df_disp['Kernel'].duplicated(), 'Kernel'] = ""


mem_df = df.groupby(["Sample Size", "Kernel"])["Memory (MiB)"].mean().reset_index(name='Memory (MiB)')
mem_df_std = df.groupby(["Sample Size", "Kernel"])["Memory (MiB)"].std().reset_index(name='std')
mem_df_std["std"] = mem_df_std["std"]/np.sqrt(10)
mem_df["std"] = mem_df_std["std"]
mem_df.fillna(0, inplace=True)



# print(time_df)
# print(mem_df)
#%%
f, ax = plt.subplots(figsize=(6, 3.5),dpi=200, facecolor="w", edgecolor="xkcd:grey")
plt.plot(time_df[time_df["Kernel"] == "SKAT_QUAD"]["Sample Size"]/1000, time_df[time_df["Kernel"] == "SKAT_QUAD"]["Time (mins)"], label="SKAT_QUAD", marker='o')
plt.plot(time_df[time_df["Kernel"] == "DIRECT"]["Sample Size"]/1000, time_df[time_df["Kernel"] == "DIRECT"]["Time (mins)"], label="QuadKAST", marker='^')
#plt.xlabel("Sample Size (k)", size=14)
plt.ylabel("Time (mins)", size=14)
#plt.title("Time (M = 100)", size=25)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 13},markerscale=2.,loc='upper right')
plt.subplots_adjust(bottom=0.17)
#plt.subplots_adjust(left=0.15)
plt.savefig("plots/time_compr_m100.png", dpi=200)
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

new_df['Time (mins)'] = new_df['Time (hrs)'] * 60
# Group by 'Sample Size' and collect 'Time (hrs)' into lists
grouped = new_df.groupby('Sample Size')['Time (mins)'].apply(list).reset_index()

# Expand lists into separate columns (Trial 1 to Trial 10)
expanded = pd.DataFrame(grouped['Time (mins)'].tolist(), index=grouped['Sample Size'])
expanded.columns = [f'Trial {i+1}' for i in range(expanded.shape[1])]
#%%
expanded.to_csv("runtimes_sample.csv")
# %%
## SNP TEST ##
# %%
import os
import pandas as pd
#%%
# Define the directory containing the data
data_directory = '50000n'

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
#%%
results_df["avg"] = results_df.iloc[:,1:].mean(axis=1)
results_df["std"] = results_df.iloc[:,1:].std(axis=1)
# %%
import matplotlib.pyplot as plt
import numpy as np
#%%
f, ax = plt.subplots(figsize=(6, 3.5))
ax.tick_params(axis='x', which='major')

xfine = np.linspace(10, 300)
tmp = curve_fit(quad_line, results_df['Window Size'], results_df['avg'])[0]

plt.errorbar(results_df['Window Size'][:-1], results_df['avg'][:-1], yerr=results_df["std"][:-1]/np.sqrt(10), linestyle='-', marker='o', markersize=7, label="QuadKAST")
# plt.plot(xfine, quad_line(xfine, tmp[0], tmp[1], tmp[2],  tmp[3],  tmp[4],), '--', linewidth=2)


# plt.xticks(ticks=[0, 1, 2, 3, 4, 5], labels=results_df["Window Size"])

plt.ylabel("Time (min)", size = 14)
plt.xlabel("Window Size", size = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=13,markerscale=2.,loc='upper left')
plt.savefig("runtimes_snp.png", dpi=200, bbox_inches='tight')
# %%
