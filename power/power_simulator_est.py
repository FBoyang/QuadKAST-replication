from sklearn.kernel_approximation import RBFSampler
from sklearn.impute import SimpleImputer
import numpy as np
import traceback
import pandas as pd
import sys
import gc
from scipy.linalg import svd
import time
from numpy.linalg import inv
from bed_reader import open_bed
import scipy
from scipy.linalg import pinvh
import fastlmmclib.quadform as qf
from chi2comb import chi2comb_cdf, ChiSquared
from sklearn.linear_model import LogisticRegression
import scipy
from numpy.core.umath_tests import inner1d
from sklearn.preprocessing import PolynomialFeatures
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
import pickle
import argparse
# from estimators_beta import *
# from utils_beta import *
sys.path.append("/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/FastKAST/")
from fastmle_res_jax import *



def simQuad(X, h2):
    X = SimpleImputer().fit_transform(X)
    N = X.shape[0]
    mapping = PolynomialFeatures((2, 2),interaction_only=True,include_bias=False)
    scaler=StandardScaler()
    X = scaler.fit_transform(X)
    Z = mapping.fit_transform(X)
    scaler=StandardScaler()
    Z = scaler.fit_transform(Z)
    D = Z.shape[1]
    eps = np.random.randn(N) * np.sqrt(1-h2)
    beta = np.random.randn(Z.shape[1]) * np.sqrt(h2)*1.0/np.sqrt(D)
    y = Z.dot(beta) + eps
    return y, Z*1.0/np.sqrt(D)


    
    

parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')

parser.add_argument('-n', '--sample', type=int, default=5000)           
parser.add_argument('-m', '--feature', type=int, default=30)      
parser.add_argument('-g', '--sigma_g', type = float, default=0.1)
parser.add_argument('-s', '--savePath', type = str, default='results/demo.pkl')
parser.add_argument('-t', '--task', type = str)
parser.add_argument('--sw', type=int, default=0)
parser.add_argument('--bfile', type=str, default='/u/project/sgss/UKBB/data/cal/filter4')
parser.add_argument('--tindex', type=int, default=1)



overwrite=False

args = parser.parse_args()
CR = 50
task = args.task
tindex = args.tindex
np.random.seed(tindex)

# set parameters
M = args.feature
N = args.sample
sigma2_g = args.sigma_g
sigma2_e = 1 - sigma2_g


superWindow = args.sw
# read genotype data
savePath = args.savePath
bfile = args.bfile
bed = bfile+'.bed'
fam = bfile+'.fam'
bim = bfile+'.bim'
gene_annot = pd.read_csv('rand_100_genes.txt')
gene_index_annot = np.loadtxt('rand_100_indices.txt')
# G = read_plink1_bin(bed, bim, fam, verbose=False)
G = open_bed(bed)
print('Finish lazy loading the genotype matrix')

bimfile = pd.read_csv(bim,delim_whitespace=True,header=None)
bimfile.columns = ['chr', 'chrpos', 'MAF', 'pos','MAJ','MIN']
famfile = pd.read_csv(fam,delim_whitespace=True,header=None)
columns = ['FID','IID','Fa','Mo','Sex','Phi']
famfile.columns = columns
# in total 22 indices, represent 22 chromosome
Windows = [] 
Posits = bimfile.iloc[:,3].values

if task == "vary_genes":
    all_results = []
    for i in range(CR*(tindex-1),CR*tindex):
        gene_row = gene_annot.iloc[i]
        gene = gene_row.gene
        CHR = gene_row.CHR
        pStart = gene_row.pStart
        pEnd = gene_row.pEnd
        Indices_bool = (bimfile.chr.astype(str)==str(CHR))&(bimfile.pos<=pEnd)&(bimfile.pos>=pStart)
        print(f'sum of Indices_bool: {np.sum(Indices_bool)}')
        iStart = np.where(Indices_bool)[0][0]
        iEnd = np.where(Indices_bool)[0][-1]
        X = G.read(index=np.s_[:,iStart:(iEnd+1)])
        X = SimpleImputer().fit_transform(X)
        y, Z = simQuad(X,sigma2_g)
        results = getfullComponentPerm(None,Z,y,VarCompEst=True)
        all_results.append(results)
        
elif task=="vary_N":
    all_results = []
    for i in range(CR*(tindex-1),CR*tindex):
        sampling = np.sort(np.random.choice(G.shape[0],N,replace=False))
        X = G.read(index=np.s_[sampling,0:M])
        # X = G[sampling,0:M]
        y, Z = simQuad(X,sigma2_g)
        results = getfullComponentPerm(None,Z,y,VarCompEst=True)
        all_results.append(results)
        
elif task=="vary_M":
    all_results = []
    for i in range(CR*(tindex-1),CR*tindex):
        iStart = int(gene_index_annot[i]) # Starting index
        print(iStart)
        X = G.read(index=np.s_[:,iStart:(iStart+M)])
        # X = G[:,iStart:(iStart+M)]
        y, Z = simQuad(X,sigma2_g)
        results = getfullComponentPerm(None,Z,y,VarCompEst=True)
        all_results.append(results)
            
elif task=="vary_h2":
    all_results = []
    for i in range(CR*(tindex-1),CR*tindex):
        iStart = int(gene_index_annot[i]) # Starting index
        X = G.read(index=np.s_[:,iStart:(iStart+M)])
        # X = G[:,iStart:(iStart+M)]
        y, Z = simQuad(X,sigma2_g)
        results = getfullComponentPerm(None,Z,y,VarCompEst=True)
        all_results.append(results)

else:
    print(f'Task {task} is undefined')   

    
with open(f'{savePath}_{tindex}.pkl', 'wb') as handle:
    pickle.dump(all_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
        
        


