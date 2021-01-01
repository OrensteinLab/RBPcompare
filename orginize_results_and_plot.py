
#TODO need to add a section to plot the results of our Z and the addition of structure information
'''
output results
# parts that in the title there is "OR" are parts that convert the txt files that
  the scripts outputs to csv files

'''
import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon as wil
from math import ceil
#check if two columns are significsnt
def IsSignif(DF):
    Pval = {}
    for i in range(len(DF.columns)):
        for j in range(i+1,len(DF.columns)):
            Pval[DF.columns[i]+'_'+DF.columns[j]] = {'X1':i,'X2':j,'p':wil(DF.iloc[:,i],DF.iloc[:,j])[1]}
    return Pval
#function to sort by the distance between significante columns
def SortByDist_sig_only(Ps):
    dis = []
    keys = []
    for key,val in Ps.items():
        if val['p']<0.05:
            keys.append(key)
            dis.append(abs(val['X1']-val['X2']))
    SortedKeys = [x for _,x in sorted(zip(dis,keys))]
    dis.sort()
    
    return SortedKeys,dis
def SortByDist(Ps):
    dis = []
    keys = []
    for key,val in Ps.items():
        
        keys.append(key)
        dis.append(abs(val['X1']-val['X2']))
    SortedKeys = [x for _,x in sorted(zip(dis,keys))]
    dis.sort()
    
    return SortedKeys,dis    
#function to plot significance        
def PltAst(DF,ax):
    M = DF.max().max()
    Ps = IsSignif(DF)
    SortedKeys, SortedDis = SortByDist(Ps)
    add = 0.05
    i=1
    index = 0
    for key in SortedKeys:
        X1 = Ps[key]['X1']
        X2 = Ps[key]['X2']
        plt.plot([X1,X1,X2,X2],[M+i*add,M+(i+1)*add,M+(i+1)*add,M+i*add],color='k')
        plt.text((X1+X2)/2,M+(i)*add,'*',horizontalalignment='center')
        i+=3

def Get_sig_level(p):
    levels = [0.05,0.01,0.001,0.0001]
    if p > levels[0]:
        return 'ns'
    if p < levels[-1]:
        return '*'*len(levels)
    for i in range(len(levels)-1):
        if p<levels[i] and p>levels[i+1]:
            return '*'*(i+1)

def PltAst1(DF,ax,winner):
    '''
    

    Parameters
    ----------
    DF : dataframe
        containes the samples that were ploted on the ax.
    ax : matplotlib axes object
        axes that the boxplots are ploted on.
    winner : string
        name of the column to plot the significance level between it to all others.

    Returns
    -------
    None.

    '''
    M = DF.max().max()
    Ps = IsSignif(DF)
    SortedKeys, SortedDis = SortByDist(Ps)
    SortedKeys2 = []
    SortedDis2 = []
    for i,key in enumerate(SortedKeys):
        if winner in key:
            SortedKeys2.append(key)
            #SortedDis2.append(SortedDis[i])
    add = 0.012
    i=1
    index = 0
    for key in SortedKeys2:
        X1 = Ps[key]['X1']
        X2 = Ps[key]['X2']
        siglevel = Get_sig_level(Ps[key]['p'])
        plt.plot([X1,X1,X2,X2],[M+i*add,M+(i+1)*add,M+(i+1)*add,M+i*add],color='k')
        if siglevel == 'ns':
            plt.text((X1+X2)/2,M+(i+1.5)*add,siglevel,horizontalalignment='center')
        else:
            plt.text((X1+X2)/2,M+(i+0.5)*add,siglevel,horizontalalignment='center')
        i+=3

#%% Figure 2

Data = pd.read_csv ('results/results_cmp_bns.csv',index_col=[0])
Data = Data.set_index([Data.index,Data.protein])
Pnames = Data['protein']
Data = Data.drop(['protein'],axis=1)
#remove the next line to use all the data available
Data = Data[['6_Z_1300_mean','6_Z_1300_max','6_R_1300_mean','5_Z_1300_mean','4_Z_1300_mean','6_Z_320_mean','6_Z_80_mean','6_Z_20_mean','6_Z_5_mean']]
Data.dropna(inplace = True)
Data.to_csv('raw_data/Figure2.csv')
used_CMP = Data.index
#average and maximum
DF = Data[['6_Z_1300_mean','6_Z_1300_max']]
DF.columns = ['mean','max']

#save for SI
DF.dropna(inplace = True)

DF.to_csv('raw_data/cmp_bns_mean_max.csv')

fig,axs = plt.subplots(nrows = 1,ncols=4,figsize = (20,5))
Pval = wil(DF['mean'],DF['max'])[1]
Min = -0.12
DF.plot.scatter('mean','max',ax = axs[0])
axs[0].plot([-0.5,0.7],[-0.5,0.7],'--',color = [1 ,0 ,0])
axs[0].set_xlim([Min-0.05,0.7])
axs[0].set_ylim([Min-0.05,0.7])
axs[0].set(title ='Mean vs Max score',
           xlabel = 'Mean',
           ylabel = 'max')
axs[0].text(0.1,0.6,'n = {0}\n p-value = {1}'.format(len(DF),ceil(Pval*1000)/1000),ha='center')

#R/Z
DF = Data[['6_Z_1300_mean','6_R_1300_mean']]            
DF.columns = ['Z','R']
Pval = wil(DF['Z'],DF['R'])[1]

#save for SI
DF.dropna(inplace = True)

DF.to_csv('raw_data/cmp_bns_R_Z.csv')

DF.plot.scatter('Z','R',ax = axs[1])
Min =-0.12
axs[1].plot([-0.5,0.7],[-0.5,0.7],'--',color = [1 ,0 ,0])
axs[1].set_xlim([Min-0.05,0.7])
axs[1].set_ylim([Min-0.05,0.7])
axs[1].set(title ='R vs Z score',
           xlabel = 'Z',
           ylabel = 'R')
axs[1].text(0.1,0.6,'n = {0}\n p-value = {1}'.format(len(DF),ceil(Pval*1000)/1000),ha='center')


#K length
DF = pd.DataFrame(Data['4_Z_1300_mean'])
DF.columns = ['K = 4']
DF['K = 5'] = Data['5_Z_1300_mean']
DF['K = 6'] = Data['6_Z_1300_mean']

#save for SI
DF.dropna(inplace = True)

DF.to_csv('raw_data/cmp_bns_K.csv')

sns.boxplot( data=DF,width = 0.4,palette = 'Set2',ax = axs[2])
sns.swarmplot(data=DF, color="crimson",ax = axs[2])

axs[2].set(title = 'differnet K',
            ylabel = 'Pearson correlation')

#concentration

C1300 = pd.DataFrame(Data['6_Z_1300_mean'])
C1300.columns = ['C = 1300']
C320 = pd.DataFrame(Data['6_Z_320_mean'])
C320.columns = ['C = 320']
C80 = pd.DataFrame(Data['6_Z_80_mean'])
C80.columns = ['C = 80']
C20 = pd.DataFrame(Data['6_Z_20_mean'])
C20.columns = ['C = 20']
C5 = pd.DataFrame(Data['6_Z_5_mean'])
C5.columns = ['C = 5']

DF = pd.concat([C5,C20,C80,C320,C1300],axis=1)
#save for SI
DF.dropna(inplace = True)
DF.to_csv('raw_data/cmp_bns_concentration.csv')

#check for statistical significance: 
Ps = IsSignif(DF)
ax = sns.boxplot( data=DF,width = 0.4,palette = 'Set2',ax = axs[3])
ax = sns.swarmplot(data=DF, color="crimson",ax = axs[3])

axs[3].set(title = 'differnet concentrations',
            ylabel = 'Pearson correlation')

plt.savefig('figures/RNAcompete_RBnS_comparsion.pdf',bbox_inches = 'tight',dpi = 350)
plt.savefig('figures/RNAcompete_RBnS_comparsion.png',bbox_inches = 'tight')
plt.show()

#%% Figure 3
#This section is for plot the comparison between the different parameters eCLIP and RBNS

Data = pd.read_csv ('results/results_clip_bns_full.csv',index_col=[0])
Data = Data.set_index([Data.index,Data.protein])
Pnames = Data['protein']
Data = Data.drop(['protein'],axis=1)
Data = Data[['4_Z_320_mean','4_Z_320_max','4_R_320_mean','5_Z_320_mean','6_Z_320_mean','4_Z_5_mean','4_Z_20_mean','4_Z_80_mean','4_Z_1300_mean']]
Data.dropna(inplace = True)
Data.to_csv('raw_data/Figure3.csv')
Used_ENC = Data.index #for later
fig,axs = plt.subplots(nrows = 1,ncols=4,figsize = (20,5))

#mean/max
DF = Data[['4_Z_320_mean','4_Z_320_max']]
DF.columns = ['mean','max']
DF.dropna(inplace = True)
#save for SI
DF.to_csv('raw_data/clip_bns/clip_bns_mean_max.csv')

Pval = wil(DF['mean'],DF['max'])[1]

Min = -0.12
DF.plot.scatter('mean','max',ax = axs[0])
axs[0].plot([0,1],[0,1],'--',color = [1 ,0 ,0])
axs[0].set_xlim([0,1])
axs[0].set_ylim([0,1])
axs[0].set(title = 'Mean vs Max score',
           xlabel = 'Mean',
           ylabel = 'Max')

axs[0].text(0.2,0.6,'n = {0}\n p-value = {1}'.format(len(DF),ceil(Pval*1000)/1000),ha='center')

#R/Z
DF = Data[['4_Z_320_mean','4_R_320_mean']]            
DF.columns = ['Z','R']
DF.dropna(inplace = True)
Pval = wil(DF['Z'],DF['R'])[1]

#save for SI
DF.to_csv('raw_data/clip_bns/clip_bns_R_Z.csv')

DF.plot.scatter('Z','R',ax = axs[1])
Min =-0.12
axs[1].plot([0,1],[0,1],'--',color = [1 ,0 ,0])
axs[1].set_xlim([0,1])
axs[1].set_ylim([0,1])
axs[1].set(title = 'R vs Z score',
           xlabel = 'Z',
           ylabel = 'R')

axs[1].text(0.2,0.6,'n = {0}\n p-value = {1}'.format(len(DF),ceil(Pval*1000)/1000),ha='center')


#K value
DF = pd.DataFrame(Data['4_Z_320_mean'])
DF.columns = ['K = 4']
DF['K = 5'] = Data['5_Z_320_mean']
DF['K = 6'] = Data['6_Z_320_mean']
DF.dropna(inplace = True)
#save for SI
DF.to_csv('raw_data/clip_bns/clip_bns_K.csv')

ax = sns.boxplot( data=DF,width = 0.4,palette = 'Set2',ax = axs[2])
ax = sns.swarmplot(data=DF, color="crimson",ax = axs[2])
# PltAst(DF,ax)
axs[2].set(title = 'clip_bns_differnet K',
           ylabel = 'AUC')

#concentration
C1300 = pd.DataFrame(Data['4_Z_1300_mean'])
C1300.columns = ['C = 1300']
C320 = pd.DataFrame(Data['4_Z_320_mean'])
C320.columns = ['C = 320']
C80 = pd.DataFrame(Data['4_Z_80_mean'])
C80.columns = ['C = 80']
C20 = pd.DataFrame(Data['4_Z_20_mean'])
C20.columns = ['C = 20']
C5 = pd.DataFrame(Data['4_Z_5_mean'])
C5.columns = ['C = 5']

DF = pd.concat([C5,C20,C80,C320,C1300],axis=1)
#check for statistical significance:
DF.dropna(inplace = True)
#save for SI
DF.to_csv('raw_data/clip_bns/clip_bns_concentration.csv')   
Ps = IsSignif(DF)
ax = sns.boxplot( data=DF,width = 0.4,palette = 'Set2',ax = axs[3])
ax = sns.swarmplot(data=DF, color="crimson",ax = axs[3])

axs[3].set(title = 'clip_bns_differnet concentrations',
           ylabel = 'AUC')

plt.savefig('figures/clip_bns_all_comparisons.png',dpi=300)      
plt.savefig('figures/clip_bns_all_comparisons.pdf',dpi=350)  
plt.show()

#%% Figure 4a
# BoxPlot for compete-compete comparison
Data = pd.read_csv ('results/results_cmp_bns.csv',index_col=[0])
RBPs = Data.protein
Data2 = pd.read_csv('results/results_cmp_cmp.csv',index_col=[0])       

Cols = list(Data2.columns[1::])

RBPs = Data2.protein
Data2.drop(['protein'],axis = 1,inplace = True)

Data2 = Data2[['z_mean','z_max','e_mean','e_max']]
Data2['best RNA Bind-n-Seq'] = Data['6_Z_1300_mean']

Data2.dropna(inplace = True)
ax = sns.boxplot( data=Data2,width = 0.4,palette = 'Set2')
ax = sns.swarmplot(data=Data2, color="crimson")
PltAst1(Data2,ax,'z_mean')

plt.title('RNA compete RNA compete',fontsize = 18)
plt.ylabel('Pearson correlation',fontsize = 18)
plt.xticks(fontsize = 18)
plt.savefig('figures/CMP-CMP.png',dpi=300)
plt.savefig('figures/CMP-CMP.pdf',dpi=350)
plt.show()
Data2['RBP'] = RBPs
Data2.to_csv('raw_data/Figure4_A.csv')
#%% FIgure 4b
# visualize the cmpt eclip and RBNS eCLIP results
Data = pd.read_csv('results/results_cmp_eclip.csv',index_col = [0])
Indcmp = [x.split('_')[1] for x in list(Data.index)]
Indclip = [x.split('_')[0] for x in list(Data.index)]
Data.index = Indclip


Data1 = pd.read_csv ('results/results_clip_bns_full.csv',index_col=[0])
# Data1 = Data1['4_Z_320_mean']
Data['BnS' ] = Data1['4_Z_320_mean']
Data['RBP'] = Data1.protein
Data['RNAompete'] = Indcmp
Data.dropna(inplace = True)
Data.to_csv('raw_data/Figure4_B.csv')

cols = ['z_mean','z_max','e_mean','e_max','BnS']
Data = Data[cols]

Ps = IsSignif(Data)

ax = sns.boxplot( data=Data[cols],width = 0.4,palette = 'Set2')
ax = sns.swarmplot(data=Data[cols], color="crimson")
PltAst1(Data,ax,'z_mean')
plt.title('RNCMPT RBnS and eCLIP comparison',fontsize=18)
plt.xticks(fontsize=15)
plt.xlabel('score type',fontsize=18)
plt.ylabel('AUC',fontsize=18)
plt.savefig('figures/RNCMPT_BnS and eCLIP comparison.png',dpi=350,bbox_inches='tight')  
plt.savefig('figures/RNCMPT_BnS and eCLIP comparison.pdf',dpi=350,bbox_inches='tight')  
plt.show()


##############################################################################
#%% The nex part will be to check how well our Z score is compared with RBNS Z score eCLIP data
Data1 = pd.read_csv ('results/results_clip_bns_full.csv',index_col=[0])
Data1 = Data1['4_Z_320_mean']

Data2 = pd.read_csv('results/Results_E_PUP_M_4.csv',index_col = [0])
IsSignif(Data2[['New_Z','Struct']])

Data2['Old_Z'] = Data1
sum(Data2.Struct>Data2.New_Z)
Data2.dropna(inplace = True)
Data2.drop(['RBP'],axis = 1,inplace = True)
temp = Data2.loc[Used_ENC]
temp.mean()
IsSignif(temp)
Data2.mean()
IsSignif(Data2)

#%% The nex part will be to check how well our Z score is compared with RBNS Z score RNAcompete data
Data = pd.read_csv ('results/results_cmp_bns.csv',index_col=[0])
Data = Data['6_Z_1300_mean']
Data2 = pd.read_csv('results/Results_E_PUP_M_6.csv',index_col = [0])
Data2.drop(['RBP'],axis = 1,inplace = True)
Data2.Struct = Data2.Struct.apply(lambda A:  float(A[1:-1].replace(',','').split()[0]))
Data2.New_Z = Data2.New_Z.apply(lambda A:  float(A[1:-1].replace(',','').split()[0]))
Data2['Old_Z'] = Data
#use only the RBPs that we used in figure 2
Data2 = Data2.loc[used_CMP]

IsSignif(Data2)
#%% BoxPlot for compete-compete comparison  OLD fig 4a

# Data2 = pd.read_csv('results/results_cmp_cmp.csv',index_col=[0])       
# # best_BNS = Data['6_Z_1300_mean']

# # Data2['BnS'] = best_BNS
# Cols = list(Data2.columns[1::])
# # Cols.sort()
# Data2.drop(['protein'],axis = 1,inplace = True)
# # Ps = IsSignif(Data2[Cols])
# Data2 = Data2[['z_mean','z_max','e_mean','e_max']]
# ax = sns.boxplot( data=Data2,width = 0.4,palette = 'Set2')
# ax = sns.swarmplot(data=Data2, color="crimson")
# PltAst1(Data2,ax,'z_mean')

# plt.title('RNA compete RNA compete',fontsize = 18)
# plt.ylabel('Pearson correlation',fontsize = 18)
# plt.xticks(fontsize = 18)
# plt.savefig('CMP-CMP.png',dpi=300)
# plt.savefig('CMP-CMP.pdf',dpi=350)
# plt.show()
        