import numpy as np
import pandas as pd
import glob
import scipy.stats
from scipy.stats import pearsonr
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




def trait_coverter(df):
    if 'CHR' in df.columns:
        df.CHR = df.CHR.astype(float).astype(int)
    if 'SNP' in df.columns:
        df.SNP = df.SNP.astype(float).astype(int)
    if 'Index' in df.columns:
        df['Index'] = df['Index'].astype(float).astype(int)
    
    annot = pd.read_csv('/home/boyang1995/research/traits/trait_info_cat.csv')
    
    for trait in np.unique(df.Trait):
        mathes = annot[annot.FILENAME==trait].TRAIT.values
        if len(mathes)!=0:
            df['Trait'][df['Trait']==trait] = mathes[0]
        else:
            df['Trait'][df['Trait']==trait] = trait

    df = df.sort_values(by=['Trait']).reset_index(drop=True)
    return df


def scatter_corr_plot(xval,yval,xlabel,ylabel,title,xval_se=None,yval_se=None,xaxis=10,yaxis=22,thresh=True,diag=False,xlim=None,ylim=None,colorMap=None,category=None,text=True):
    figure(figsize=(8, 6), dpi=200)
    corr, _ = pearsonr((xval),(yval))
    corr = round(corr, 3)
    # plt.style.use('ggplot')
    # plt.grid(False)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)
    plt.title(title, fontsize=30)
    if text:
        plt.text(xaxis,yaxis, f'$\\rho = $ {corr}',fontsize=20)
    logy=yval
    logx=xval
    if thresh:
        intersect_index = (logy<=-np.log10(5e-8))
        # plt.scatter(-np.log10(p_all_whole[('array','P')]+1e-160), -np.log10(p_all_whole[('pc40 array','P')]+1e-160), alpha=0.5, c='red')
        if (xval_se is not None) or (yval_se is not None):
            plt.errorbar(logx[intersect_index], logy[intersect_index],xerr=xval_se[intersect_index],yerr=yval_se[intersect_index],alpha=0.5,c='blue',fmt="o")
            plt.errorbar(logx[~intersect_index], logy[~intersect_index],xerr=xval_se[~intersect_index],yerr=yval_se[~intersect_index],alpha=0.5,c='red',fmt="o")
        else:
            plt.scatter(logx[intersect_index], logy[intersect_index], alpha=0.5, c='blue')
            plt.scatter(logx[(~intersect_index)], logy[(~intersect_index)], alpha=0.5, c='red')
        plt.axhline(y=-np.log10(5e-8), color='r', linestyle='--')
        plt.axvline(x=-np.log10(5e-8), color='r', linestyle='--')
    else:
        print(logx, logy)
        if (xval_se is not None) or (yval_se is not None):
            if colorMap is None:
                plt.errorbar(logx,logy,xerr=xval_se,yerr=yval_se,alpha=0.5,c='blue',fmt="o")
            else:
                print(colorMap)
                scatter_colors = np.array([colorMap[cat] for cat in category])
                for i in range(len(logx)):
                    plt.errorbar(logx[i],logy[i],xerr=xval_se[i],yerr=yval_se[i],alpha=0.5,color=colorMap[category[i]],fmt="o")
                
        else:
            if colorMap is None:
                plt.scatter(logx, logy, alpha=0.5, c='blue')
            else:
                scatter_colors = np.array([colorMap[cat] for cat in category])
                plt.scatter(logx, logy, alpha=0.5, c=scatter_colors)
    if colorMap is not None:
        for cat in np.unique(category):
            plt.scatter(logx[category==cat], logy[category==cat], alpha=0.5, c=scatter_colors[category==cat],label=cat)
        plt.legend(fontsize=15)
    if diag: # add diagonal line
        diag_array=np.linspace(-10,max(np.max(logx),np.max(logy))*1.05,100)
        plt.plot(diag_array,diag_array,color='r',linestyle='--')
    if xlim is not None:
        plt.xlim(xlim[0],xlim[1])
    if ylim is not None:
        plt.ylim(ylim[0],ylim[1])
    plt.locator_params(axis='x', nbins=6)
    plt.locator_params(axis='y', nbins=6)
    plt.show()