'''
this is a script for comparing RNAcompete results to RNA Bind-n-Seq results
'''

import pandas as pd
import numpy as np

import scipy.stats as stat
import pickle
from utiles import *
dircmp = '/data/yaronore/RNAcompete/'
dirbns = '/data/yaronore/RNA_bind_n_seq/'
#dircmp=''
#dirbns=''
'''
#convert the scores table to dictionary for faster activation
#==============================================================================
def ScoresToDic(Scores,exp,proteins,ez):
    #exp should be 'cmpt' for RNAcompete 7-mers scores or 'bns' for ....
    if exp == 'cmpt':
        
        Scores.index = Scores['7mer']
        
        #delete the columns of the sets we don't need
        Cnames =list( Scores.columns)
        Cnames1=[]
        for name in Cnames:
            if 'setB' in name or 'setAB' in name:
                Cnames1.append(name)

        Scores = Scores.drop(['7mer']+Cnames1,axis=1)
        #delete the proteins colomns we don't need
        cols2keep=[]
        for p in proteins:
            cols2keep.append(p+'_'+ez+'_setA')
        Scores = Scores[cols2keep]
        #add the missing 7mers
        df2 = pd.DataFrame(Scores.mean()).T
        df2.index = ['GAAGAGC']
        Scores = Scores.append(df2)
        df3 = pd.DataFrame(Scores.mean()).T
        df3.index = ['GCTCTTC']
        Scores = Scores.append(df3)        
        return Scores.to_dict()
#    else:
        
#calculate the appropriate value for a single sequence    
#==============================================================================        
def Seq2Val(Seq,Max_Mean,SDic,k):
#    print(len(SDic))
    Seq = Seq.replace('U','T')
    Max_Mean = Max_Mean.upper()
    tmpScores = []
    for Start in range(len(Seq)-k+1):
        try:
            tmpScores.append(SDic[Seq[Start:Start+k]])
        except:
            tmpScores.append(SDic['mean'])    
    return [max(tmpScores) , np.mean(tmpScores)]
#    if Max_Mean =='MAX':
#        return max(tmpScores)
#    else:
#        return np.mean(tmpScores)

def ReadScores(protein,RZ,k):
    # dirbns = '/data/yaronore/RNA_bind_n_seq/overlap/'
    dirbns = '/data/eitamar/TechComp/tsv_files/'
    Scores = pd.read_table(dirbns+p+'_'+RZ+'_{0}.tsv'.format(k))
    Scores.index = Scores.iloc[:,0]
    L = Scores.columns
    Scores = Scores.drop([L[0]],axis=1)
    df2 = pd.DataFrame(Scores.mean(axis=0)).T
    df2.index=['mean']
    Scores = Scores.append(df2)    
    return Scores
'''
###############################################################################    
####################################  MAIN  ###################################
###############################################################################

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
        
        
        
K = list(range(4,8))
RZ=['R','Z']


results={}
for pcode,p in PDic.items():
    print(p)
    results[pcode]={}
    f = open(dircmp+pcode+'.txt.sequences_B.RNAcontext.clamp')
    org = []
    Seqs = []
    for line in f:
        spl = line.split(' ')
        org.append(float(spl[0]))
        Seqs.append(spl[1])
  
    
    for k in K:
        results[pcode]['{0}'.format(k)]={}
        for rz in RZ:
            results[pcode]['{0}'.format(k)][rz]={}
            try:
                Scores = ReadScores(p,rz,k)
            except:
                continue
            SDic = Scores.to_dict()
            concs = list(Scores)   
            for conc in concs:
                results[pcode]['{0}'.format(k)][rz][conc]={}
                new_max = []
                new_mean = []
                for seq in Seqs:
                    tmp = Seq2Val(seq,SDic[conc],k)
                    new_max.append(tmp[0])
                    new_mean.append(tmp[1])
                
                results[pcode]['{0}'.format(k)][rz][conc]={'max':stat.pearsonr(org,new_max),'mean': stat.pearsonr(org,new_mean)}  

                with open('results/RNAcmpt_RNAbns_comparison.txt','wb') as fb:    
                    pickle.dump(results,fb)
                
Pro=[]
Pro2=[]
Ks=[]
RZs=[]
Concs=[]
mM=[]



dfs={}
for P_key,P_val in results.items():
    Cols=[]
    scores=[]
    for K_key,K_val in P_val.items():
        for RZ_key,RZ_val in K_val.items():
            for C_key,C_val in RZ_val.items():
                if 'nM' in C_key:
                    Cols.append(K_key+'_'+RZ_key+'_'+C_key[0:-3]+'_'+'max')
                    Cols.append(K_key+'_'+RZ_key+'_'+C_key[0:-3]+'_'+'mean')
                else:                    
                    Cols.append(K_key+'_'+RZ_key+'_'+C_key+'_'+'max')
                    Cols.append(K_key+'_'+RZ_key+'_'+C_key+'_'+'mean')
                scores.append(C_val['max'][0])
                scores.append(C_val['mean'][0])
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
DF.to_csv('results/results_cmp_bns.csv')

        
        