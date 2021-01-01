'''
this is a script for comaring RNAcompete results to RNAcompete results
'''

import pandas as pd
import numpy as np
#from utiles import GetScore,IntToBase,CharToint, BreakToKmers
#import joblib
import scipy.stats as stat
import pickle
from utiles import *
dircmp = '/data/yaronore/RNAcompete/'
# dircmp=''

#convert the scores table to dictionary for faster activation
#==============================================================================
# def ScoresToDic(Scores,exp,proteins,ez):
#     #exp should be 'cmpt' for RNAcompete 7-mers scores or 'bns' for ....
#     if exp == 'cmpt':
        
#         Scores.index = Scores['7mer']
        
#         #delete the columns of the sets we don't need
#         Cnames =list( Scores.columns)
#         Cnames1=[]
#         for name in Cnames:
#             if 'setB' in name or 'setAB' in name:
#                 Cnames1.append(name)

#         Scores = Scores.drop(['7mer']+Cnames1,axis=1)
#         #delete the proteins colomns we don't need
#         cols2keep=[]
#         for p in proteins:
#             cols2keep.append(p+'_'+ez+'_setA')
#         Scores = Scores[cols2keep]
#         #add the missing 7mers
#         df2 = pd.DataFrame(Scores.mean()).T
#         df2.index = ['GAAGAGC']
#         Scores = Scores.append(df2)
#         df3 = pd.DataFrame(Scores.mean()).T
#         df3.index = ['GCTCTTC']
#         Scores = Scores.append(df3)        
#         return Scores.to_dict()
#    else:
        
#calculate the appropriate value for a single sequence    
#==============================================================================        
# def Seq2Val(Seq,Max_Mean,SDic,k):
#     Max_Mean = Max_Mean.upper()
#     tmpScores = []
#     for Start in range(len(Seq)-k):
#         tmpScores.append(SDic[Seq[Start:Start+k]])
#     return [max(tmpScores) , np.mean(tmpScores)]

    
#%% ===========================================================================

Codes = pd.read_table('BnS_CMPT.txt',sep=' ')
#create dictionary for protein's names
PDic={}
L=list(Codes.index)
i=0
for cmp in Codes['CMP']:
    splited = cmp.split(',')
    for code in splited:
        PDic[code] = L[i]
    i+=1
    
proteins=[]
for cmp in Codes['CMP']:
    splited = cmp.split(',')
    for i in splited:
        proteins.append(i)
ZScores = pd.read_table(dircmp+'z_scores.txt')
ZDic = ScoresToDic(ZScores,'cmpt',proteins,'z')
del ZScores
EScores = pd.read_table(dircmp+'e_scores.txt')
EDic = ScoresToDic(EScores,'cmpt',proteins,'e')
del EScores

results={}
for p in proteins:
    # p='RNCMPT00001'
    f = open(dircmp+p+'.txt.sequences_B.RNAcontext.clamp')
    org = []
    new_max_z = []
    new_mean_z = []
    new_max_e = []
    new_mean_e = []
    for line in f:
        
        spl = line.split(' ')
        org.append(float(spl[0]))
        #zscore
        tmp = Seq2Val(spl[1],ZDic[p+'_z_setA'],7)
        new_max_z.append(tmp[0])
        new_mean_z.append(tmp[1])
        #escore
        tmp = Seq2Val(spl[1],EDic[p+'_e_setA'],7)
        new_max_e.append(tmp[0])
        new_mean_e.append(tmp[1])
        
    pears_max = stat.pearsonr(org,new_max_z)
    pears_mean = stat.pearsonr(org,new_mean_z)  
    
    results[p]={'z_max':stat.pearsonr(org,new_max_z),
           'z_mean': stat.pearsonr(org,new_mean_z),
           'e_mean': stat.pearsonr(org,new_mean_e),
           'e_max': stat.pearsonr(org,new_max_e)}    
        
    with open('results/RNAcmpt_RNAcmpt comparison.txt','wb') as fb:    
        pickle.dump(results,fb)


dfs={}
for P_key,P_val in results.items():
    Cols=[]
    scores=[]
    for key,val in P_val.items():
        Cols.append(key)
        scores.append(val[0])

    DF = pd.DataFrame(scores).T
    DF.columns = Cols
    DF.index = [P_key]
    DF['protein'] = PDic[P_key]
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
DF.to_csv('results/results_cmp_cmp.csv')




        
        
