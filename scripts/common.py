###
 # @author Wankun Deng
 # @email dengwankun@hotmail.com
 # @create date 2024-07-05 09:31:58
 # @modify date 2024-07-05 09:31:58
 # @desc [description]
###
import warnings
warnings.filterwarnings('ignore')

import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
from multiprocessing import Pool

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import scipy
from statsmodels.stats.multitest import multipletests
from scipy.stats import zscore

import scanpy as sc
import anndata
import scvelo as scv

from joblib import dump, load

from sklearn.metrics import *
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier



os.chdir('/home/wdeng3/workspace/Codespace/scRTE/scripts/')

dataset='../data/all_datasets/AD_HS_00001.cell_exp.txt'
rmsk_f='../../universal_data/rmsk/rmsk_GRCh38.txt'
dataset_umap='../data/all_datasets/AD_HS_00001.cell_umap.txt'
gtf='../../universal_data/ref/GRCh38/gencode.v43.basic.annotation.gtf'
input_path='../data/all_datasets/'

gene_type={}
for line in open(gtf):
    if not line.startswith('#'):
        info=line.strip().split('\t')
        if info[2]=='gene':
            split_='gene_name "' if 'gene_name' in info[-1] else 'gene_id "'
            gene_name=info[-1].split(split_)[1].split('"')[0]
            if 'lncRNA' in info[-1] :
                gene_type[gene_name]='lncRNA'
            elif 'protein_coding' in info[-1]:
                gene_type[gene_name]='protein_coding'
            elif 'pseudogene' in info[-1]:
                gene_type[gene_name]='pseudogene'
            else:
                gene_type[gene_name]='Others'
            
rmsk=pd.read_csv(rmsk_f,sep='\t')
rmsk['repFamily']=[x.replace('?','') for x in rmsk['repFamily']]
rmsk['repLen']=rmsk['genoEnd']-rmsk['genoStart']
classification=rmsk.loc[rmsk['repClass'].isin(['SINE', 'LINE', 'LTR']),['repName','repClass','repFamily']].drop_duplicates(ignore_index=True)
tmp=rmsk.iloc[:,[10,11]].drop_duplicates()
te_cls=dict(zip(tmp['repName'],tmp['repClass']))
tmp=rmsk.iloc[:,[10,11,12]].drop_duplicates()
tmp=tmp.loc[tmp['repClass'].isin(['SINE','LINE','LTR']),:]
te_fam=dict(zip(tmp['repName'],tmp['repFamily']))
fam_te=defaultdict(list)
for x in te_fam:
    fam_te[te_fam[x]].append(x)
gene_type.update(te_cls)
rtes=rmsk['repName'].unique()

genes=list(gene_type.keys())
genes_rep=[x.replace('_','.').replace('-','.') for x in genes]
file_list=[x for x in os.listdir(input_path) if x.endswith('.cell_exp.txt')]

def get_dataset(dataset):
    print(f'loading data: {dataset} \n')
    dt_ls=[x for x in file_list if x.startswith(dataset)]
    print('Reading %s \n'%os.path.join(input_path,dt_ls[0]))
    cell_exp=pd.read_table(os.path.join(input_path,dt_ls[0]),index_col=0)
    if len(dt_ls) >1:
        cell_umap=pd.read_table(f'{input_path}/{dataset}.1.cell_umap.txt',index_col=0)
        for i in range(1,len(dt_ls)):
            print('Reading %s \n'%os.path.join(input_path,dt_ls[i]))
            cell_exp=pd.concat([cell_exp,pd.read_table(f'{input_path}/'+dt_ls[i],index_col=0)])
            cell_umap=pd.concat([cell_umap,pd.read_table(f'{input_path}/{dataset}.{i}.cell_umap.txt',index_col=0)])
    else:
        cell_umap=pd.read_table(f'{input_path}/'+dataset+'.cell_umap.txt',index_col=0)
    
    cell_umap['predicted.celltype'] = cell_umap['predicted.celltype'].replace(
        'Opc', 'OPC')
    for i in range(cell_umap.shape[0]):
        if cell_umap.iloc[i,1] =='Stage_0':
            cell_umap.iloc[i,1]='Control'
        if cell_umap.iloc[i,1] !='Control' and  not cell_umap.iloc[i,1].startswith('Stage'):
            cell_umap.iloc[i,1]=cell_umap.iloc[i,7].split('_')[0]
    
    colnames=cell_exp.columns
    repl_colnames=[]
    for x in colnames:
        if '.' not in x or x not in genes_rep:
            repl_colnames.append(x)
        else:
            repl_colnames.append(genes[genes_rep.index(x)])
    cell_exp.columns=repl_colnames
    print(f'Done loading: {dataset} \n')
    return [cell_exp,dataset,cell_umap]

def get_dataset_nomerge(dataset):
    print(f'loading data: {dataset} \n')
    cell_exp=pd.read_table(f'{input_path}/{dataset}.cell_exp.txt',index_col=0)
    cell_umap=pd.read_table(f'{input_path}/{dataset}.cell_umap.txt',index_col=0)
    
    cell_umap['predicted.celltype'] = cell_umap['predicted.celltype'].replace(
        'Opc', 'OPC')
    for i in range(cell_umap.shape[0]):
        if cell_umap.iloc[i,1] =='Stage_0':
            cell_umap.iloc[i,1]='Control'
        if cell_umap.iloc[i,1] !='Control' and  not cell_umap.iloc[i,1].startswith('Stage'):
            cell_umap.iloc[i,1]=cell_umap.iloc[i,7].split('_')[0]
    
    colnames=cell_exp.columns
    repl_colnames=[]
    for x in colnames:
        if '.' not in x or x not in genes_rep:
            repl_colnames.append(x)
        else:
            repl_colnames.append(genes[genes_rep.index(x)])
    cell_exp.columns=repl_colnames
    print(f'Done loading: {dataset} \n')
    return [cell_exp,dataset,cell_umap]

pool=Pool(15)
datasets=[x.replace('.cell_exp.txt','') for x in file_list]
results=pool.map(get_dataset_nomerge,datasets)                                                                                                                                                                    
pool.close()
pool.join()
all_dfs={}
all_cell_umaps={}
def load_sc_data():
    for ret in results: 
        df,dataset,cell_umap=ret
        cell_umaps=df.iloc[:,-2:]
        df=np.expm1(df.iloc[:,:-2])
        df=pd.concat([df,cell_umaps],axis=1)
        all_dfs[dataset]=df.copy()
        all_cell_umaps[dataset]=cell_umap.copy()

    results=None