"""
this script should be run after the csv files of the matricis created
to run this proparly you should use the command:
    python prot_All.py --bns 1 --pup 0 -e 1
bns - should the script process the RNA bind_n_seq sequences and structure
    files and create the scores of the Kmers?
pup - should the analysis be based on the entire structural information or
     paired/un-paired information (5 by K OR 2 by K structural context matrix)
enrich - should the script use enrichment of the structure of every Kmer or
    mean matrix? 
"""

import pandas as pd
import optparse
import os
from multiprocessing import Pool
import pickle
import numpy as np
#==============================================================================

def getOptions():
    parser = optparse.OptionParser()    

    parser.add_option("--bns", action="store", type="int", default=1,
                      help = "process RBnS files? yes = 1, No = 0")  

    parser.add_option("--pup", action="store", type="int", default=1,
                      help = "paired / un-paired matrix or full? paired un -paired = 1, Full = 0")
    
    parser.add_option("-e", "--enrich", action="store", type="int", default=1,
                      help = "create enrichment matirx or mean matrix? mean = 0, enrichment = 1")     
    
    parser.add_option("-c", "--collectRes", action="store", type="int", default=0,
                      help = "collect all results to one txt file? yes = 1 no = 0")  
    
    parser.add_option("-m", "--mean", action="store", type="int", default=1,
                      help = "use mean over all nucleotides? 0 = no , otherwise = yes")

    parser.add_option("--eCLIP", action="store", type="int", default = 1,
                      help = "predict on eCLIP files? 0 = no , otherwise = yes")
    
    parser.add_option("--cmp", action="store", type="int", default=0,
                      help = "predict on RNAcompete files? 0 = no , otherwise = yes")
    
    parser.add_option("-K", action="store", type="int", default=4,
                      help = "length of kmers")   
    
    (options, args) = parser.parse_args()    
    return options  

#==============================================================================
#define function so could run n parrallel
def Run_Command(cmd):
        os.system(cmd)
        return('DONE '+cmd)

#==============================================================================        

###############################################################################    
####################################  MAIN  ###################################
###############################################################################


dircmp = '/data/yaronore/RNAcompete/'
dirbns = '/data/yaronore/RNA_bind_n_seq/'
direclip = '/data/yaronore/eCLIP/'
direclip = '/data/eitamar/eCLIP_files/'

options = getOptions()
PDic={}
#create all possible kmers
os.system(f'./CreateKmers {options.K} > {options.K}.txt')

if options.cmp:
    #load paires of RNA bind n seq and RNAcompete
    Codes = pd.read_table('BnS_CMPT.txt',sep=' ')
    L=list(Codes.index)
    i=0
    for cmp in Codes['CMP']:
        splited = cmp.split(',')
        for code in splited:
            PDic[code] = L[i]
        i+=1
        
#load the eCLIP paires with RNA Bind-n-Seq
if options.eCLIP:    
    Codes = pd.read_table('eCLIP_BnS.txt',sep = '\t') 
    for i in range(len(Codes)):
        PDic[Codes['eCLIP'].loc[i]] = Codes['protein'].loc[i]

#create a list of commands that need to be executed          
CMD_create_scores = []     
CMD_give_scores = []
DONE=[]
for cmp,bns in PDic.items():
    if (bns not in DONE) and options.bns:
        #os.system("python count_prot.py {0}".format(bns))
        CMD_create_scores.append("python CreateScoresFiles.py {0} --pup {1} -e {2} -m {3} -K {4}"
                   .format(bns,options.pup,options.enrich,options.mean,options.K))
        DONE.append(bns)
       
    CMD_give_scores.append("python GiveScores.py {0} {1} --pup {2} -e {3} -m {4} -K {5}"
               .format(bns,cmp,options.pup,options.enrich,options.mean,options.K))

        
#so that the Create Scores files will run first    
#CMD.sort()

#run up to "parallel_exp_num" parrallel commands
if options.bns:
    parallel_exp_num = 10  # 7
    with Pool(parallel_exp_num) as p:  # 7 processes at a time
        reslist = [p.apply_async(Run_Command, (CMD_create_scores[i],)) for i in range(len(CMD_create_scores))]
        for result in reslist:
            #pass
            print(result.get())
 
parallel_exp_num = 20  # 7
with Pool(parallel_exp_num) as p:  # 7 processes at a time
    reslist = [p.apply_async(Run_Command, (CMD_give_scores[i],)) for i in range(len(CMD_give_scores))]
    for result in reslist:
        #pass
        print(result.get())
        
#collect all results to one file
if options.collectRes :
    #try to open if a file exists
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
        
    Resfile = 'Results_{0}_{1}_{2}.txt'.format(E,pup,M)
    try:
        with open(Resfile,'rb') as fb:
            Res = pickle.load(fb)
    except:
        Res = {}
        
    for cmp,bns in PDic.items():
        try: #because not all RBPs have all concentrations
            with open('{0}_{1}_{2}_{3}_{4}.txt'.format(cmp,E,pup,M,options.K),'rb') as fb:
                Res[cmp] = pickle.load(fb)
            os.system('rm {0}_{1}_{2}_{3}_{4}.txt'.format(cmp,E,pup,M,options.K))
        except:
            continue
    Resfile = 'Results_{0}_{1}_{2}_{3}.txt'.format(E,pup,M,options.K)  
    with open(Resfile,'wb') as fb:
        pickle.dump(Res,fb)



    # pup = 'PUP'
    # K = 4
    # E='E'
    # M = 'M'
    # conc = 320
    # with open('results/Results_{0}_{1}_{2}_{3}.txt'.format(E,pup,M,options.K),'rb') as fb:
    #     Res = pickle.load(fb)
    
    #orginize the ressults to CSV
    if options.eCLIP:
        conc = 320
    else:
        conc = 1300
    cods = []
    RBPs = []
    struct = []
    our_z = []
    for code,res in Res.items():
        cods.append(code)
        our_z.append(res[conc]['Z'])
        struct.append(res[conc]['dot_mean'])
        RBPs.append(res['BnS'])
    DF = pd.DataFrame(RBPs,index= cods)        
    DF.columns = ['RBP']
    DF['New_Z'] = our_z
    DF['Struct'] = struct
    DF.to_csv('results/Results_{0}_{1}_{2}_{3}.csv'.format(E,pup,M,options.K))









    
    
    
    
