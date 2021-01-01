'''
This is a script uses the k-mers scores as computed by RNA Bind-n-Seq developers 
to predict in-vivo binding in eCLIP experiments

To use this script you must have:
* A folder containing all .tsv files (that contain the k-mers scores)
* Folder containing all eCLIP sequences files
'''

import pandas as pd
import numpy as np
from sklearn import metrics
import scipy.stats as stat
import pickle
from utiles import ReadScores,Score_eclip,Seq2Val
dircmp = '/data/yaronore/RNAcompete/'
direclip = '/data/eitamar/TechComp/eCLIP_files/'
dirbns = '/data/eitamar/TechComp/tsv_files/'


#%% ===========================================================================
#reads the overlap between the tecnologies 
Codes = pd.read_table('eCLIP_BnS.txt')
proteins = list(Codes['protein'])
#create dictionary for protein's names
PDic={}
temp = list(Codes['eCLIP'])
for i in range(len(temp)):
    PDic[temp[i]] = proteins[i]




#define some parameters      
K = list(range(4,8))
RZ=['R','Z']

results={}
for clip_code,p in PDic.items():
    print(p)
    results[clip_code]={}

    for k in K:
        results[clip_code]['{0}'.format(k)]={}
        for rz in RZ:
            results[clip_code]['{0}'.format(k)][rz]={}
            
            #read RBNS file
            try:
                Scores = ReadScores(p,rz,k)
                print('loaded scores')
            except:
                #if there is no scores file for the specific protein or K
                continue
            SDic = Scores.to_dict()
            concs = list(Scores)
            #end of reading RBNS file
            
            for conc in concs:
#                results[clip_code]['{0}'.format(k)][rz][conc]={}

                # load negatives and score them
                Neg_max,Neg_mean = Score_eclip(clip_code+'.bed.control.ext.hg19.fa',SDic[conc],k) 
            
                # load positives and score them
                Pos_max,Pos_mean = Score_eclip(clip_code+'.bed.ext.hg19.fa',SDic[conc],k)  
                
                #compute AUC for maximum aggregation function
                truelabels = [1]*len(Pos_max)+[0]*len(Neg_max) 
                fpr,tpr,th = metrics.roc_curve(truelabels,Pos_max+Neg_max,pos_label=1)
                tmpMAX = metrics.auc(fpr,tpr)
                ROCmax = {'tpr':tpr.copy(),'fpr':fpr.copy(),'thresh':th.copy()}
                #compute AUC for average aggregation function
                fpr,tpr,th = metrics.roc_curve(truelabels,Pos_mean+Neg_mean,pos_label=1)
                tmpMEAN = metrics.auc(fpr,tpr)               
                ROCmean = {'tpr':tpr.copy(),'fpr':fpr.copy(),'thresh':th.copy()}
                #store the results
                results[clip_code]['{0}'.format(k)][rz][conc]={'mean':tmpMEAN,'max':tmpMAX ,'ROC':{'mean':ROCmean,'max':ROCmax}}
                
                #save the results in a dictionary file
                with open('results/eCLIP_RNAbns_comparison_new.txt','wb') as fb:    
                    pickle.dump(results,fb)
                

#convert the results to a csv file and saves it
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
                scores.append(C_val['max'])
                scores.append(C_val['mean'])
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
DF.to_csv('results/results_clip_bns_full.csv')

 