'''
utiles file
'''

import numpy as np
import pandas as pd

def Seq2Val(Seq,SDic,k):
#    print(len(SDic))
    Seq = Seq.replace('U','T')
#    Max_Mean = Max_Mean.upper()
    tmpScores = []
    for Start in range(len(Seq)-k+1):
        try:
            tmpScores.append(SDic[Seq[Start:Start+k]])
        except:
            if '\n' in Seq[Start:Start+k]:
                break
            tmpScores.append(SDic['mean'])  
    if len(tmpScores) == 0:
        tmpScores.append(np.nan)  
    return [max(tmpScores) , np.mean(tmpScores)]


def Score_eclip(clip_Fname,Kmers_Scores,k):
    direclip = '/data/yaronore/eCLIP/'
    direclip = '/data/eitamar/TechComp/eCLIP_files/'
    # direclip =''
    f = open(direclip+clip_Fname)
    Seq_score_mean = []
    Seq_score_max = []
#    print(Kmers_Scores['mean'])
    for line in f:
        if '>' in line:
            continue
        line = line[150:-150].upper()
        Seq_score=Seq2Val(line,Kmers_Scores,k)
        if np.nan in Seq_score:
            continue
        Seq_score_max.append(Seq_score[0])
        Seq_score_mean.append(Seq_score[1])
    return Seq_score_max,Seq_score_mean

#reads the RNA Bind-n-Seq scores files
#=============================================================================
def ReadScores(protein,RZ,k):
    print(RZ+'_{0}.tsv'.format(k))
    # dirbns = '/data/yaronore/RNA_bind_n_seq/overlap/'
    dirbns = '/data/eitamar/TechComp/tsv_files/'

    print(dirbns+protein+'_'+RZ+'_{0}.tsv'.format(k))
    Scores = pd.read_table(dirbns+protein+'_'+RZ+'_{0}.tsv'.format(k))
    Scores.index = Scores.iloc[:,0]
    L = Scores.columns
    Scores = Scores.drop([L[0]],axis=1)
    df2 = pd.DataFrame(Scores.mean(axis=0)).T
    df2.index=['mean']
    Scores = Scores.append(df2)    
    return Scores

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
        df4 = pd.DataFrame(Scores.mean()).T
        df4.index = ['mean']
        Scores = Scores.append(df4) 
        return Scores.to_dict()