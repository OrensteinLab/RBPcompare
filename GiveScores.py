"""
This script will use the scores files created using "CreateScoresFiles_v1.py",
score the "cmp" and compute the pearson correlation between the scores computed
by this script and the original RNAcompete scores

if the "cmp" argumant is eCLIP code for protein the calculation will be of AUC
and not pearson correlation


                = = = = = = = arguments = = = = = = =
                
pup - should the analysis be based on the entire structural information or
     paired/un-paired information (5 by K OR 2 by K structural context matrix)
enrich - should the script use enrichment of the structure of every Kmer or
    mean matrix? 
    

"""

import scipy.stats as stat
import pickle
import numpy as np
import sys
import math
import optparse
from sklearn import metrics
from utiles import *
missing=0
#RNA bind n seq protein name
protein = sys.argv[1] 
# protein = 'TARDBP'
# cmp = 'ENCFF869RDN'  
#RNAcompete \ eCLIP experiment identifier
cmp = sys.argv[2]

dir_eCLIP    = '/data/yaronore/eCLIP/'
dir_eCLIP    = '/data/eitamar/TechComp/eCLIP_files/'
dir_CMPT     = '/data/yaronore/RNAcompete/'
s_dir = 'scores_files/'


#==============================================================================
def getOptions():
    parser = optparse.OptionParser()    

    parser.add_option("--pup", action="store", type="int", default=0,
                      help = "paired / un-paired matrix or full? paired un -paired = 1, Full = 0")
    
    parser.add_option("-e", "--enrich", action="store", type="int", default=1,
                      help = "create enrichment matirx or mean matrix? mean = 0, enrichment = 1") 

    parser.add_option("-K", action="store", type="int", default=6,
                      help = "The Kmers length")    
    
    parser.add_option("-m", "--mean", action="store", type="int", default=1,
                      help = "use mean over all nucleotides? 0 = no , otherwise = yes")
    
    parser.add_option("-c", "--concentration", action="store", type="int", default=1300,
                      help = "Whihch RBnS RBP concentration k-mers scores to use?")

    (options, args) = parser.parse_args()    
    return options 


#==============================================================================
    
def ReadCompScores(fname):
    f = open(fname,'r')
    org=[]
    for line in f:
        org.append(float(line.split(' ')[0]))
    return org

#==============================================================================
    
def mean2(x):
    y = np.sum(x) / np.size(x);
    return y

#==============================================================================
    
def corr2(a,b):
    a = a - mean2(a)
    b = b - mean2(b)

    r = (a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum());
    return r

#==============================================================================
    
def Seq2ValZ(seq,Dic):
    tmpcorr = []
    for i in range(len(seq)-K+1):
        try:
            tmpcorr.append(Dic[seq[i:i+K]]['Z']) 
        except:
            tmpcorr.append(Dic['mean']['Z'])    
    if len(tmpcorr) == 0:
        tmpcorr.append(np.nan)        
    return np.mean(tmpcorr)

#==============================================================================

def Seq2Val1(seq,Mat,Dic):
    tmpcorr = []
    Mat1 = np.asarray(Mat).astype(float)
    for i in range(len(seq)-K+1):
        tmpcorr_per=[]
        try:
            mat2score = Dic[seq[i:i+K]]['Mat']
        except:
            continue
            mat2score = np.asarray([[0.2]*K]*5)    
        section_mat = Mat1[:,i:i+K]
        # tmpcorr.append(np.sum(np.multiply(section_mat,mat2score)))
        for j in range(K):
            tmpcorr_per.append(stat.pearsonr(mat2score[:,j],section_mat[:,j])[0])
        Corr = np.mean(tmpcorr_per)
        tmpcorr.append(Corr*Dic[seq[i:i+K]]['Z'])
    return np.mean(tmpcorr)

#==============================================================================
    
# def Seq2Val(seq,Mat,Dic):
#     tmpcorr = []
#     Mat1 = np.asarray(Mat).astype(float)
#     for i in range(len(seq)-K+1):
#         try:
#             mat2score = Dic[seq[i:i+K]]['Mat']
#         except:
#             continue
#             mat2score = np.asarray([[0.2]*K]*5)    
#         section_mat = Mat1[:,i:i+K]
        
#         tmp = section_mat[:2,:].copy()
#         tmp[1,:] = sum(section_mat[1:,])
#         # tmpcorr.append(np.sum(np.multiply(section_mat,mat2score)))
#         if options.pup:
#             tmpcorr.append(Dic[seq[i:i+K]]['Z']*corr2(tmp,mat2score))
#         else:
#             tmpcorr.append(Dic[seq[i:i+K]]['Z']*corr2(section_mat,mat2score)) 
#     return np.mean(tmpcorr)

#==============================================================================
def Seq2Val2(seq,Mat,Dic):
    tmpcorr = []
    Mat1 = np.asarray(Mat).astype(float)
    for i in range(len(seq)-K+1):
        try:
            mat2score = Dic[seq[i:i+K]]['Mat']
        except:
            continue
            mat2score = np.asarray([[0.2]*K]*5)    
        section_mat = Mat1[:,i:i+K]
        

        # tmpcorr.append(np.sum(np.multiply(section_mat,mat2score)))
        if options.pup:
            tmp = section_mat[:2,:].copy()
            tmp[1,:] = sum(section_mat[1:,])
            tmpcorr.append(Dic[seq[i:i+K]]['Z']*((tmp-mat2score)**2).mean(axis=None))
        else:
            tmpcorr.append(Dic[seq[i:i+K]]['Z']*((section_mat-mat2score)**2).mean(axis=None))
            
    return np.mean(tmpcorr)

#==============================================================================
#that is the function to use for now    
def Seq2ValM(seq,Mat,Dic):
    tmpcorr = []
    Mat1 = np.asarray(Mat).astype(float)
    for i in range(len(seq)-K+1):
        try:
            mat2score = Dic[seq[i:i+K]]['Mat']
        except:
            #this is for k mers with N (unidentified nucleotide)
            mat2score = Dic['mean']['Mat']    
        
        section_mat = Mat1[:,i:i+K]

        # tmpcorr.append(np.sum(np.multiply(section_mat,mat2score)))
        if options.pup:
            tmp = section_mat[:2,:].copy()
            tmp[1,:] = sum(section_mat[1:,])
            tmpcorr.append(np.sum(np.multiply(tmp.mean(axis=1),mat2score),axis=None))
            # tmpcorr.append(corr2(tmp,mat2score))
        else:
            tmpcorr.append(np.sum(np.multiply(section_mat.mean(axis=1),mat2score),axis=None))
        if len(tmpcorr) == 0:
            tmpcorr.append(np.nan)
        #tmpcorr = tmpcorr.sort(rever)    
    return np.mean(tmpcorr)

#==============================================================================
    
# def Seq2ValM_pre(seq,Mat,Dic):
#     precent = 0.6
#     tmpcorr = []
#     Mat1 = np.asarray(Mat).astype(float)
#     for i in range(len(seq)-K+1):
#         try:
#             mat2score = Dic[seq[i:i+K]]['Mat']
#         except:
#             continue
#             mat2score = np.asarray([[0.2]*K]*5)    
#         section_mat = Mat1[:,i:i+K]

#         # tmpcorr.append(np.sum(np.multiply(section_mat,mat2score)))
#         if options.pup:
#             tmp = section_mat[:2,:].copy()
#             tmp[1,:] = sum(section_mat[1:,])
#             tmpcorr.append(np.sum(np.multiply(tmp.mean(axis=1),mat2score),axis=None))
#             # tmpcorr.append(corr2(tmp,mat2score))
#         else:
#             tmpcorr.append(np.sum(np.multiply(section_mat.mean(axis=1),mat2score),axis=None))
        
#         tmpcorr.sort(reverse=True)
#         if int(precent*len(tmpcorr))>0:
#             tmpcorr = tmpcorr[:int(precent*len(tmpcorr))]
        
        
#     return np.mean(tmpcorr)
#==============================================================================
#this function is fot options.e==0
#the idea is a dot product between mean probabilities (0,1) and multiply it by Z    
def Seq2Val_Ztdot(seq,Mat,Dic):
    tmpcorr = []
    Mat1 = np.asarray(Mat).astype(float)
    for i in range(len(seq)-K+1):
        try:
            mat2score = Dic[seq[i:i+K]]['Mat']
        except:
            continue
            mat2score = np.asarray([[0.2]*K]*5)    
        section_mat = Mat1[:,i:i+K]

        # tmpcorr.append(np.sum(np.multiply(section_mat,mat2score)))
        if options.pup:
            tmp = section_mat[:2,:].copy()
            tmp[1,:] = sum(section_mat[1:,])
            tmpcorr.append(Dic[seq[i:i+K]]['Z']*np.sum(np.multiply(tmp.mean(axis=1),mat2score),axis=None))
            # tmpcorr.append(corr2(tmp,mat2score))
        else:
            tmpcorr.append(Dic[seq[i:i+K]]['Z']*np.sum(np.multiply(section_mat.mean(axis=1),mat2score),axis=None))
        
        #tmpcorr = tmpcorr.sort(rever)    
    return np.mean(tmpcorr)

#==============================================================================
    
def eCLIP2Scores(code,Dic,positive):
    #load the files
    if positive:
        sfile = open(dir_eCLIP+'{0}.bed.ext.hg19.fa'.format(code),'r') 
        pfile = open(dir_eCLIP+'{0}.bed.ext.hg19.fa.combined_profile.txt.sho.txt'.format(code),'r') 
    else:
        sfile = open(dir_eCLIP+'{0}.bed.control.ext.hg19.fa'.format(code),'r') 
        pfile = open(dir_eCLIP+'{0}.bed.control.ext.hg19.fa.combined_profile.txt.sho.txt'.format(code),'r') 

        
        
    scores_s=[]
    scores_z=[]
    sindex = 0
    #go through the sequnces
    for sline in sfile:
        if '>' in sline:
            sindex+=1
            continue

        pindex = 0
        probs = []
        #read the probabilities
        for pline in pfile:
            if pindex==0:
                pindex+=1
                continue
            probs.append(pline.split())
            probs[pindex-1] = [float(i) for i in probs[pindex-1]]
            pindex+=1
            if pindex == 6:
                break

        #assign a score to the sequence
        probs = np.asarray(probs)
        # print(probs)
        # print(sline)
        probs = probs[:,150:-150]
        sline = sline[150:-150]
        sline = sline.upper()
        # print(sline)
        temp_s = Seq2ValM(sline,probs,Dic)
        temp_z = Seq2ValZ(sline,Dic)
        # temp_z = Seq2Val(sline,Dic,K)[1]
        if np.isnan(temp_z):
            continue
        scores_s.append(temp_s)
        scores_z.append(temp_z)

    return {'Z':scores_z,'S':scores_s}


     
###############################################################################    
####################################  MAIN  ###################################
###############################################################################


#define the parameters for the program
options = getOptions()
K  = options.K

if options.enrich:
    E = 'E'
else:
    E = 'noE'
    
if options.pup:
    pup = 'PUP'
else:
    pup = 'Full'
if options.mean:
    M = 'M'
else:
    M = 'noM'
    
Fname = s_dir + '{0}_Scores_{1}_{2}_{3}_{4}.txt'.format(protein,E,pup,M,K)

#first, check if it is RNAcompete or eCLIP files:
if 'ENC' in cmp:
    eCLIP = True
else:
    eCLIP = False
    #read the original RNAcompete scores
    org = ReadCompScores(dir_CMPT+'{0}.txt.sequences_B.RNAcontext.clamp'.format(cmp))
    [float(i) for i in org]

#load the scores file (matrices)
with open(Fname,'rb') as fb:
    Scores = pickle.load(fb)       

#check for the highest concentrarion and use it only:
conces=[]
for key,val in Scores.items():
    if key == 'input':
        continue
    conces.append(int(key))

if not eCLIP:
    conc = max(conces)
    conc = 1300
else:
    conc = 320
# conc = options.concentration
print('{0}_{1}_{2}'.format(conc,protein,cmp))
          

Pears={}              
Dic = Scores['{0}'.format(conc)]


#open the file with the sequences
#if RNAcompete
if not eCLIP:
    seqf = open(dir_CMPT+'{0}.txt.annotations_B.RNAcontext'.format(cmp))

    index= 0
    Mat=[]
    seqscore = []
    score2 = []
    scoreZ = []
    for line in seqf:
        if index%6 == 0:
            seq = line[1:-1].replace('U','T')
            index +=1
            continue
        else:
            Mat.append(np.asarray(line[1:-1].split('\t')).transpose())
            
            
        if index%6==5:
            #first type of score calculation
            # if options.mean:
            #     score = Seq2Val_Ztdot(seq,Mat,Dic)
            probs = np.asarray(Mat)
            score = Seq2ValM(seq,probs,Dic)
            seqscore.append(score)
            #calculation 
            scoreZ.append(Seq2ValZ(seq,Dic))
            Mat=[]
        index+=1
    
    Pears[conc] = {} 
    Pears[conc]['dot_mean'] = stat.pearsonr(org,seqscore)
    # Pears[conc]['crr_mean_struct'] = stat.pearsonr(org,score2)
    #save the next two for checking no errors were made
    Pears[conc]['Z'] = stat.pearsonr(org,scoreZ)
    Pears['BnS'] = protein
    

# the next part is for calculating the eCLIP scores
if eCLIP:
    #positive part
    pos = eCLIP2Scores(cmp,Dic,positive = True)
    
   
    #negative part
    neg = eCLIP2Scores(cmp,Dic,positive = False)
    
    #calculate AUC - Z score
    truelabels = [1]*len(pos['Z'])+[0]*len(neg['Z'])
    fpr,tpr,th = metrics.roc_curve(truelabels,pos['Z']+neg['Z'],pos_label=1)
    AUC_Z = metrics.auc(fpr,tpr)
    
    #calculate AUC struct
    truelabels = [1]*len(pos['S'])+[0]*len(neg['S'])
    fpr,tpr,th = metrics.roc_curve(truelabels,pos['S']+neg['S'],pos_label=1)
    AUC_S = metrics.auc(fpr,tpr)
    
    #for saving the results
    Pears[conc] = {} 
    Pears[conc]['dot_mean'] = AUC_S
    #save the next two for checking no errors were made
    Pears[conc]['Z'] = AUC_Z
    Pears['BnS'] = protein
    Pears['ROC'] = {'tpr':tpr,'fpr':fpr,'thresh':th}
           
    
with open('{0}_{1}_{2}_{3}_{4}.txt'.format(cmp,E,pup,M,K),'wb') as fb:
   pickle.dump(Pears, fb)   
    






