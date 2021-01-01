"""
this script is aimed to create a file with a dictionary with the scores of all
Kmers
The matirecs can represent paired / un-paired probabilities or enrichmrent 
(dependt on -e value) 
pup - should the analysis be based on the entire structural information or
     paired/un-paired information (5 by K OR 2 by K structural context matrix)
enrich - should the script use enrichment of the structure of every Kmer or
    mean matrix? 
"""
import pickle
import pandas as pd
import numpy as np
import sys
import optparse

s_dir = 'scores_files/'

cb1 = 1
#==============================================================================
def getOptions():
    parser = optparse.OptionParser()    

    parser.add_option("--pup", action="store", type="int", default=1,
                      help = "paired / un-paired matrix or full? paired un -paired = 1, Full = 0")
    
    parser.add_option("-e", "--enrich", action="store", type="int", default=1,
                      help = "create enrichment matirx or mean matrix? mean = 0, enrichment = 1")   
    
    parser.add_option("-K", action="store", type="int", default=4,
                  help = "The Kmers length") 

    parser.add_option("-m", "--mean", action="store", type="int", default=1,
                      help = "use mean over all nucleotides in kmer? 0 = no , otherwise = yes")
    
    (options, args) = parser.parse_args()    
    return options  
#==============================================================================
    
def CorrectP(mat,Count,K):
    mat = mat.astype(float)
    mat[0,:] = int(Count)-np.sum(mat[1:,:],axis=0)
    return mat

#==============================================================================
    
def ReadScores(protein,conc,K):
    nm = 'nM'
    if conc == 'input':
        nm = ''
    else:
        conc = conc.split(' ')[0]
    f = open(s_dir + '{0}_{1}{2}_{3}.csv'.format(protein,conc,nm,K))
    Dic={}
    flag=1
    for line in f: 
        if flag:
            flag=0
            continue
        
        line = line[:-1]
        SPL = line.split(',')
        Kint = SPL[0]
        Count = SPL[1]
        SPL = SPL[2::]
        mat = np.asarray(SPL)
        mat = mat.reshape((5,K)) 
        mat = CorrectP(mat,Count,K)
        Dic[Kmers[int(Kint)]] = {'Count':Count,'Mat':mat}
    return Dic

#==============================================================================
    
def Mat2PUP(mat):
    temp = mat[0:2,:].copy()
    temp[1,:] = sum(mat[1:,:])
    return temp
    
#==============================================================================
def Read_Kmers(k):
    Kmers = []
    f = open(f'{k}.txt','r')
    for line in f:
        Kmers.append(line[:-1])
    return Kmers
###############################################################################    
####################################  MAIN  ###################################
###############################################################################
options = getOptions()
K = options.K   

if cb1:
    protein = sys.argv[1]
    # with open('/data/eitamar/TechComp/RNAplfold_orginize/ProtAndConc1.txt','rb') as fb: #TOSO 
    #         Dic = pickle.load(fb)
    with open('ProtANDConc.txt','rb') as fb: #TOSO 
            Dic = pickle.load(fb)    
    # mer = pd.read_table('/data/yaronore/RNA_bind_n_seq/overlap/TIA1_Z_6.tsv')
    # mer = pd.read_table(f'{K}k.txt',header = None)
    
else:
    protein='TIA1'
    #get the Kmers in acending order.        
    mer = pd.read_table('TIA1_Z_6.tsv')
    with open('C:\\Users\\Eitamar\\Google Drive\\Masters\\TechCompare\\ProtAndConc1.txt','rb') as fb:
            Dic = pickle.load(fb)

     
# mers = mer.iloc[:,0]
# Kmers = list(mers.sort_values())
Kmers = Read_Kmers(K)


#prepare the dictionary for the Scores, dictionary for each concentration
Scores = {}
conc = Dic[protein]
for c in conc:
    c = c.split(' ')[0]
    if c=='0':
        continue
    try:
        Scores[c]={}
        Scores[c] = ReadScores(protein,c,K)     
    except:
        del Scores[c]
    
#change the input matrices to float so in the next stage we could divide the
#other matrices by it
for Kint,val in Scores['input'].items():
    Scores['input'][Kint]['Mat'] = Scores['input'][Kint]['Mat'].astype(float)
    if options.pup:
        Scores['input'][Kint]['Mat'] = Mat2PUP(Scores['input'][Kint]['Mat'])
    if options.mean:
        Scores['input'][Kint]['Mat'] = Scores['input'][Kint]['Mat'].mean(axis=1)
        

for c,val in Scores.items():
    if c=='input':
        continue
    for Kint,val1 in val.items():
        val1['Mat'] = val1['Mat'].astype(float)
        if options.pup:
            val1['Mat'] = Mat2PUP(val1['Mat'])
        
        if options.mean:
            #from 5 by 6 to 5 by 1
            val1['Mat'] = val1['Mat'].mean(axis=1)
            
        if options.enrich:
             #mat        =mat/mat_input
            val1['Mat'] = np.divide(val1['Mat'],Scores['input'][Kint]['Mat'])
        else:
            #mean matrix
            val1['Mat'] =val1['Mat']/int(val1['Count'])

        val1['Z'] =  int(val1['Count'])/int(Scores['input'][Kint]['Count'])  


#this part add mean Z and matrix (in case there is a k-mer missing)
for c,val in Scores.items():
    if c =='input':
        continue
    MeanZ = 0
    MeanMat = np.zeros_like(val1['Mat'])
    ind = 0
    for Kint,val2 in val.items():
        ind+=1
        MeanZ += val2['Z']
        MeanMat = MeanMat+val2['Mat']
    Scores[c].update({'mean':{'Z':MeanZ/(4**K)   ,'Mat': MeanMat/(4**K) }})

    
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

with open(Fname,'wb') as fb:
    pickle.dump(Scores,fb)      



