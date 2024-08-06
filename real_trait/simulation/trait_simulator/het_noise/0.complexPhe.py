#%%
import numpy as np
import pandas as pd
import os
import sys
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from scipy import stats
from bed_reader import open_bed
import matplotlib.pyplot as plt

# %%

np.random.seed(1)
seed = int(sys.argv[1])
print(f'seed is {seed}')

def sim_baseline(G, cauRatio, var_g, N=None, M=None,chunksize=500):
    if M is None:
        M = G.shape[1]
    if N is None:
        N = G.shape[0]
    
    Mcau = int(M*cauRatio)
    M_indices = np.sort(np.random.choice(M, Mcau, replace=False))
    additive_eff = np.random.randn(Mcau)*np.sqrt(var_g/Mcau)
    
    
    chunks = int(np.ceil(Mcau*1.0/chunksize))
    
    y_base = np.zeros(N)
    
    for c in range(chunks):
        c_idx = M_indices[c*chunksize:(c+1)*chunksize]
        Xt = G.read(np.s_[:,c_idx])
        imputer=SimpleImputer(missing_values=np.nan, strategy='mean')
        Xt=imputer.fit_transform(Xt)
        
        scaler=StandardScaler()
        Xt=scaler.fit_transform(Xt)
        
        y_base += Xt@additive_eff[c*chunksize:(c+1)*chunksize]
    
    # y_base += np.random.randn(N)*np.sqrt(var_e)
    return y_base
    

def ivrt_transform(phe_target_orig):
    phe_target = phe_target_orig.copy()
    non_nan_data = phe_target.iloc[:, 2].dropna()
    # print(non_nan_data)
    nan_data = phe_target.iloc[:, 2].isna()
    # print(f'nan_data # is: {np.sum(nan_data)}')

    # Rank and transform only non-NaN data
    ranks = stats.rankdata(non_nan_data, method='average')

    quantiles = (ranks - 0.5) / len(non_nan_data)

    normalized_non_nan_data = stats.norm.ppf(quantiles)
    # Replace the non-NaN values in the original DataFrame with the transformed values
    valid_indices = np.where(~nan_data)[0]
    phe_target.iloc[valid_indices,2 ] = normalized_non_nan_data

    return phe_target

#%%
# covarPath='/u/project/sriram/alipazok/first_40_pcs.txt'

# covar_df = pd.read_csv(covarPath,sep=' ')

np.random.seed(1)


var_g=0.3
var_hom_e=0.3
var_het_e=0.4
cauRatio=0.001
pars=[f'cau_{cauRatio}_h2_{var_g}']


gen="/u/project/sgss/UKBB/data/cal/filter4.bed"

famdf = pd.read_csv("/u/project/sgss/UKBB/data/cal/filter4.fam",sep=' ',header=None).iloc[:,:2].values



G = open_bed(gen)
# pars=['cau_0.1_h2_0.1']
rootPhePath='/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG/code/robustness/Nat.Gen/Simulation_small/'
#%%


for par in pars:
    np.random.seed(seed)
    
    # pheno_df_orig = pd.read_csv(f'{rootPhePath}{par}/{i}.pheno',delimiter=' ')
    # sim_baseline(G, cauRatio, var_g, var_e, N=None, M=None,chunksize=500)
    N = G.shape[0]
    M = G.shape[1]
    pheno_df_orig = sim_baseline(G, cauRatio=cauRatio, var_g=var_g, chunksize=500)
    # N = pheno_df_orig.shape[0]
    Indicator = np.zeros(N)
    het_samples = np.random.choice(N,int(0.5*N),replace=False)
    Indicator[het_samples] = 1
    
    eff_noise = np.sqrt(var_het_e*Indicator + var_hom_e)
    pheno_df = np.concatenate((famdf, pheno_df_orig.reshape(-1,1)),axis=1)
    pheno_df = pd.DataFrame(data=pheno_df,columns=['FID','IID','pheno'])
    # Indicator_var=(Indicator-np.mean(Indicator))/np.std(Indicator)
    pheno_df.iloc[:,2] = pheno_df.iloc[:,2] + eff_noise*np.random.randn(N)
    # pheno_df.columns = ['FID','IID','pheno']
    
    pheno_df_ivrt = ivrt_transform(pheno_df)
    # pheno_df_ivrt = ivrt_transform(pheno_df_orig)
    
    
    dir_path = f'traits/{par}'
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    # # pheno_df_ivrt.to_csv(f'{dir_path}/{i}.orig.ivrt.pheno',sep=' ',index=None)
    pheno_df_ivrt.to_csv(f'{dir_path}/{seed}.het.pheno',sep=' ',index=None)
        
        
# %%
