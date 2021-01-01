'''
this is a script for comaring RNAcompete results to RNAcompete results
'''

import pandas as pd
import numpy as np
from sklearn import metrics
import scipy.stats as stat
import pickle
from utiles import *
dircmp = '/data/yaronore/RNAcompete/'
dirbns = '/data/yaronore/RNA_bind_n_seq/'
direclip = '/data/yaronore/eCLIP/'
#dircmp=''
#dirbns=''
#direclip = ''

#convert the scores table to dictionary for faster activation
#==============================================================================
def ScoresToDic(Scores,exp,proteins,ez):
    #exp should be 'cmpt' for RNAcompete 7-mers scores or 'bns' for ....
    if exp == 'cmpt':
        
        Scores.index = Scores['7mer'].str.replace('U','T')
        
        #delete the columns of the sets we don't need
        Cnames =list( Scores.columns)
        Cnames1=[]
        for name in Cnames:
            if 'setB' in name:
                Cnames1.append(name)

        Scores = Scores.drop(['7mer']+Cnames1,axis=1)
        #delete the proteins colomns we don't need
        cols2keep=[]
        for p in proteins:
            cols2keep.append(p+'_'+ez+'_setAB')
        Scores = Scores[cols2keep]
        #add the missing 7mers
        df2 = pd.DataFrame(Scores.mean()).T
        df2.index = ['GAAGAGC']
        Scores = Scores.append(df2)
        df3 = pd.DataFrame(Scores.mean()).T
        df3.index = ['GCUCUUC']
        Scores = Scores.append(df3)  
        df4 = pd.DataFrame(Scores.mean()).T
        df4.index = ['mean']
        Scores = Scores.append(df4)   
        
        return Scores.to_dict()


###############################################################################    
####################################  MAIN  ###################################
###############################################################################

#load the codes
PList = []
proteins = []
f = open('CMPT_eCLIP.txt','r')
for line in f:
    PList.append(line[:-1])
    proteins.append(line.split('_')[1][:-1])
    
proteins = list(set(proteins))

ZScores = pd.read_table(dircmp+'z_scores.txt')
ZDic = ScoresToDic(ZScores,'cmpt',proteins,'z')
del ZScores
EScores = pd.read_table(dircmp+'e_scores.txt')
EDic = ScoresToDic(EScores,'cmpt',proteins,'e')
del EScores
    
    
K=7
results={}
i=0
for code in PList:
    i+=1
    
    print('pare {0} from {1}'.format(i,len(PList)+1))
    print(code)
    results[code]={}
    splited = code.split('_')
    clip_code = splited[0]
    CMPT_code = splited[1]
    #Zscore
    # load negatives
    Neg_max,Neg_mean = Score_eclip(clip_code+'.bed.control.ext.hg19.fa',ZDic[CMPT_code+'_z_setAB'],K) 

    # load positives
    Pos_max,Pos_mean = Score_eclip(clip_code+'.bed.ext.hg19.fa',ZDic[CMPT_code+'_z_setAB'],K)     
    
    truelabels = [1]*len(Pos_max)+[0]*len(Neg_max)
    
    fpr,tpr,th = metrics.roc_curve(truelabels,Pos_max+Neg_max,pos_label=1)
    results[code]['z_max'] = metrics.auc(fpr,tpr)
    fpr,tpr,th = metrics.roc_curve(truelabels,Pos_mean+Neg_mean,pos_label=1)
    results[code]['z_mean'] = metrics.auc(fpr,tpr)    
    
    #Escore
    # load negatives
    Neg_max,Neg_mean = Score_eclip(clip_code+'.bed.control.ext.hg19.fa',EDic[CMPT_code+'_e_setAB'],K) 

    # load positives
    Pos_max,Pos_mean = Score_eclip(clip_code+'.bed.ext.hg19.fa',EDic[CMPT_code+'_e_setAB'],K)     
    
    truelabels = [1]*len(Pos_max)+[0]*len(Neg_max)
    
    fpr,tpr,th = metrics.roc_curve(truelabels,Pos_max+Neg_max,pos_label=1)
    results[code]['e_max'] = metrics.auc(fpr,tpr)
    fpr,tpr,th = metrics.roc_curve(truelabels,Pos_mean+Neg_mean,pos_label=1)
    results[code]['e_mean'] = metrics.auc(fpr,tpr)

    

    with open('results/eCLIP_CMPT_compare.txt','wb') as fb:    
        pickle.dump(results,fb)
                

dfs={}
for P_key,P_val in results.items():
    Cols=[]
    scores=[]
    for key,val in P_val.items():
        Cols.append(key)
        scores.append(val)

    DF = pd.DataFrame(scores).T
    DF.columns = Cols
    DF.index = [P_key]
#    DF['protein'] = PDic[P_key]
    dfs[P_key] = DF

i=-1
for key,val in dfs.items():
    i+=1
    if not i:
        DF = val
    else:
        DF = pd.concat([DF,val],axis=0)
Cols = list(DF)        
Cols1 = [Cols[len(Cols)-1]]
for i in range(len(Cols)-1):
    Cols1.append(Cols[i])
DF=DF[Cols1]
Cols = list(DF.columns)
Cols.sort()
DF = DF[Cols]
DF.to_csv('results/results_cmp_eclip.csv')
        
        