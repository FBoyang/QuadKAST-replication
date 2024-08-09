from sklearn.kernel_approximation import RBFSampler
import numpy as np
import pandas as pd
import traceback
from scipy.optimize import minimize
from bed_reader import open_bed
import sys
import gc
from scipy.linalg import svd
import time
from numpy.linalg import inv
import scipy
from scipy.linalg import pinvh
import fastlmmclib.quadform as qf
from chi2comb import chi2comb_cdf, ChiSquared
from sklearn.linear_model import LogisticRegression
import scipy
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from numpy.core.umath_tests import inner1d
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
import pickle
import argparse
from scipy import stats
sys.path.append("/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/FastKAST/")
from fastmle_res_jax import *


def direct(geno_matrix_in):
    N=geno_matrix_in.shape[0]
    M=geno_matrix_in.shape[1]
    exact = np.zeros((N, int((M*(M+1))/2)))
    s = 0
    for i in range(M):
        for j in range(i, M):
            feature = geno_matrix_in[:,i]*geno_matrix_in[:,j]
            if j != i:
                feature *= np.sqrt(2)
            exact[:,s] = feature
            s += 1
    exact_standard = stats.zscore(exact)

    # remove nan columns
    col_mean = np.nanmean(exact_standard, axis=0)
    inds = np.where(np.isnan(exact_standard))
    exact_standard[inds] = np.take(col_mean, inds[1])
    #exact_standard = exact_standard[:,~np.all(np.isnan(exact_standard), axis=0)]
    return exact_standard

parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')

parser.add_argument('--bfile', type=str, default='/u/project/sgss/UKBB/data/cal/filter4')
parser.add_argument('--phenPath', default='/u/project/sriram/alipazok/Data/new_ukkb_phenotypes/', help='Phenotype file. Required')
parser.add_argument('--tindex', type=int, default=1)
parser.add_argument('--sw', type=int, default=2, help='The superwindow is set to a multiple of the set dimension at both ends, default is 2')
parser.add_argument('--savePath', type=str, default='/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/down_stream/results1/')

args = parser.parse_args()
tindex=args.tindex
bfile=args.bfile
superWindow=args.sw

### load the meta information
sig_meta = pd.read_csv('/u/scratch/a/aanand2/QuadKAST/quadkast_sig_genes.csv')
entry = sig_meta.iloc[tindex-1]
trait, gene, pStart, pEnd, Start, End, Morig, Porig = entry

savePath = args.savePath
bfile = args.bfile
bed = bfile+'.bed'
fam = bfile+'.fam'
bim = bfile+'.bim'

# G = read_plink1_bin(bed, bim, fam, verbose=False)
G = open_bed(bed)
print('Finish lazy loading the genotype matrix')

bimfile = pd.read_csv(bim,delim_whitespace=True,header=None)
bimfile.columns = ['chr', 'chrpos', 'MAF', 'pos','MAJ','MIN']
famfile = pd.read_csv(fam,delim_whitespace=True,header=None)
columns = ['FID','IID','Fa','Mo','Sex','Phi']
famfile.columns = columns
Windows = [] 


CHR = bimfile.iloc[Start].chr
## prepare covariates
covar=f'{args.phenPath}/{trait}.covar'
covarfile = pd.read_csv(covar,delim_whitespace=True)
assert covarfile.iloc[:,0].equals(famfile.FID)
covarfile = covarfile.iloc[:,2:]

phe=f'{args.phenPath}/{trait}.pheno'
Y = pd.read_csv(phe,delim_whitespace=True).iloc[:,-1].values
Yeffect = (Y!=-9)&(~np.isnan(Y))
Y = Y[Yeffect]
covarfile = covarfile[Yeffect]


start = Start
end = End+1
wlen = end - start
c = 2-G.read(index=np.s_[Yeffect,max(0,start-superWindow*wlen):min(G.shape[1],end+superWindow*wlen)])
c = np.concatenate((c,covarfile),axis=1)
x = 2-G.read(index=np.s_[Yeffect,start:end])
nanfilter=~np.isnan(c).any(axis=1)

c = c[nanfilter]
c = np.unique(c, axis=1, return_index=False)
x = x[nanfilter]
x = np.unique(x, axis=1, return_index=False)
y = Y[nanfilter]
scaler=StandardScaler()
c = scaler.fit_transform(c)
scaler=StandardScaler()
x = scaler.fit_transform(x)



### Interaction only ####
mapping = PolynomialFeatures((2, 2),interaction_only=True,include_bias=False)
Z = mapping.fit_transform(x)
# Z = direct(x)
scaler=StandardScaler()
Z = scaler.fit_transform(Z)
D = Z.shape[1]
Z = Z*1.0/np.sqrt(D)

### Get variance components
results = getfullComponentPerm(c,Z,y,VarCompEst=True,center=True)
results['pval-orig']=Porig
print(f"p-orig: {Porig} vs p-recal: {results['pval']}")
results['meta']=[trait,gene,CHR,pStart,pEnd]

### Get the feature importance
y_res = y-LinearRegression().fit(c,y).predict(c)

model = sm.OLS(y_res, Z)
model = model.fit()

# Getting coefficients
results['OLS-pvals']=model.pvalues
results['OLS-params']=model.params

with open(f'{savePath}_{trait}_{gene}.pkl', 'wb') as handle:
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

### Include the self-interaction ###

mapping = PolynomialFeatures((2, 2),interaction_only=False,include_bias=False)
Z = mapping.fit_transform(x)
# Z = direct(x)
scaler=StandardScaler()
Z = scaler.fit_transform(Z)
D = Z.shape[1]
Z = Z*1.0/np.sqrt(D)

### Get variance components
results = getfullComponentPerm(c,Z,y,VarCompEst=True,center=True)
results['pval-orig']=Porig
print(f"p-orig: {Porig} vs p-recal-inter: {results['pval']}")
results['meta']=[trait,gene,CHR,pStart,pEnd]

### Get the feature importance
y_res = y-LinearRegression().fit(c,y).predict(c)

model = sm.OLS(y_res, Z)
model = model.fit()

# Getting coefficients
results['OLS-pvals']=model.pvalues
results['OLS-params']=model.params

with open(f'{savePath}{trait}_{gene}.inter.pkl', 'wb') as handle:
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)



