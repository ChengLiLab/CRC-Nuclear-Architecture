#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:07:00 2021

@author: liyaru
"""
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import wilcoxon
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt

files=os.listdir("./out_SDOC_1")
name=[]
all_SDOC={}

for fl in files:
    if fl[-4:]==".tsv":
        name.append(fl[0:4]+"_1")
        f=open("./out_SDOC_1/"+fl)
        lines=f.readlines()
        for l in lines:
            l=l.strip().split()
            tad=l[0]+"_"+l[1]+"_"+l[2]
            if tad not in all_SDOC:
                all_SDOC[tad]=[]
            all_SDOC[tad].append(float(l[-1]))

files=os.listdir("./out_SDOC_2")
for fl in files:
    if fl[-4:]==".tsv":
        name.append(fl[0:4]+"_2")
        f=open("./out_SDOC_2/"+fl)
        lines=f.readlines()
        for l in lines:
            l=l.strip().split()
            tad=l[0]+"_"+l[1]+"_"+l[2]
            if tad not in all_SDOC:
                all_SDOC[tad]=[]
            all_SDOC[tad].append(float(l[-1]))

TAD2=[]
f=open("/home/liyaru/software/temp/SDOC_all/tad_all_boundary/all_TAD_2")
lines=f.readlines()
for l in lines:
    l=l.strip().split("\t")
    tad=l[0]+"_"+l[1]+"_"+l[2]
    TAD2.append(tad)

all_SDOC2={}
for a in all_SDOC:
    if len(all_SDOC[a])==12 and (a in TAD2):
        all_SDOC2[a]=all_SDOC[a]

df=pd.DataFrame(all_SDOC2)
df=df.T
df.columns=name
df2=df[["04CN_1","04CN_2","05CN_1","05CN_2","14CN_1","14CN_2",
       "04CT_1","04CT_2","05CT_1","05CT_2","14CT_1","14CT_2"]]

p_value=[]
mean_T_N=[]
for i in df2.index:
    N=df2.loc[i][0:6].tolist()
    T=df2.loc[i][6:].tolist()
    
    #p_value.append(stats.mannwhitneyu(N,T)[1])
    #两独立样本 曼-惠特尼U检验（Mann-Whitney U test），又称曼-惠特尼秩和检验，
    #可以看作是对两均值之差的参数检验方式的T检验或相应的大样本正态检验的代用品    
    
    p_value.append(wilcoxon(N,T)[1])
    #它适用于T检验中的成对比较，但并不要求成对数据之差di服从正态分布，只要求对称分布即可。
    #检验成对观测数据之差是否来自均值为0的总体
    
    #p_value.append(ttest_rel(N,T)[1])
    #配对T test
    
    d=np.mean(T)-np.mean(N)
    mean_T_N.append(d)
    
df2["p_value"]=p_value
df2["mean_T_N"]=mean_T_N
df2.to_csv("./analysis_data/NT_SDOC_P_value.csv")

df3=df2.loc[df2["p_value"]<0.05]
#df3=df2.loc[df2["p_value"]<0.10]

df_up=df3.loc[df3["mean_T_N"]>0]
df_down=df3.loc[df3["mean_T_N"]<0]

df_up.to_csv(path_or_buf="./analysis_data/NT_up_SDOC_TAD")
df_down.to_csv(path_or_buf="./analysis_data/NT_down_SDOC_TAD")


#看与TSS overlap
data_file1={}
data_file2={}

f=open("./analysis_data/TSS_hg19_ganjb.txt")
lines=f.readlines()
for l in lines:
    l=l.strip().split("\t")
    if l[0] not in data_file1:
        data_file1[l[0]]=[]
    data_file1[l[0]].append(l)    
    
for i in df_down.index:
#for i in df_down.index:
    i=i.split("_")
    if i[0] not in data_file2:
        data_file2[i[0]]=[]
    data_file2[i[0]].append(i[0:3])    

common_chr=sorted(list(set(data_file1)&set(data_file2)))

#得到SDOC改变的gene列表
down_SDOC_gene=[]
for chrom in common_chr:
    for l1 in data_file1[chrom]:
        s1=int(l1[1])
        e1=int(l1[2])
        for l2 in data_file2[chrom]:
            s2=int(l2[1])
            e2=int(l2[2])
            if e1>s2 and e2>s1:
                down_SDOC_gene.append(l1[5])
down_SDOC_gene2=list(set(down_SDOC_gene))               

#看这些基因在自己数据中的表达量
exp=pd.read_csv("./analysis_data/all_sample_gene_expr_normalized2.txt",sep="\t")
exp=exp.set_index(['V1'])
exp=exp[['CRC-04-N', 'CRC-05-N','CRC-14-N', 
         'CRC-04-T', 'CRC-05-T','CRC-14-T']]

'''
exp["mean_N"]=exp[['CRC-04-N', 'CRC-05-N','CRC-14-N']].mean(axis=1) 
exp["mean_T"]=exp[['CRC-04-T', 'CRC-05-T','CRC-14-T']].mean(axis=1) 
#自己的数据 整体的表达量 N和T 是没有差异的
print(ttest_rel(exp["mean_N"],exp["mean_T"]))
'''

down_SDOC_gene3=[]
for e in down_SDOC_gene2:
    if e in exp.index:
        down_SDOC_gene3.append(e)
              
diff_SDOC_exp=exp.loc[down_SDOC_gene3,:]
   
diff_SDOC_exp["mean_N"]=diff_SDOC_exp[['CRC-04-N', 'CRC-05-N','CRC-14-N']].mean(axis=1) 
diff_SDOC_exp["mean_T"]=diff_SDOC_exp[['CRC-04-T', 'CRC-05-T','CRC-14-T']].mean(axis=1) 

diff_SDOC_exp.to_csv(path_or_buf="./analysis_data/down_SDOC_gene_exp",sep="\t")
#diff_SDOC_exp.to_csv(path_or_buf="./analysis_data/up_SDOC_gene_exp",sep="\t")

print(wilcoxon(diff_SDOC_exp["mean_N"],diff_SDOC_exp["mean_T"],alternative="less"))
print(wilcoxon(diff_SDOC_exp["mean_N"],diff_SDOC_exp["mean_T"],alternative="greater"))
print(ttest_rel(diff_SDOC_exp["mean_N"],diff_SDOC_exp["mean_T"]))
plt.figure()
f=plt.boxplot([diff_SDOC_exp["mean_N"],diff_SDOC_exp["mean_T"]],labels=["N","T"],widths=0.1,positions=[0,0.15])


#看这些基因在TCGA中的表达量
df=pd.read_csv("/media/liyaru/LYR/colon/27TCGA/TCGA/tcga_normalized2.txt",sep="\t")
df['Ensembl gene ID'] = df['Ensembl_ID'].apply(lambda x: x.split(".")[0])

df2=pd.read_csv("/media/liyaru/LYR/colon/27TCGA/TCGA/geneID_transfer.txt",sep="\t")

tcga=pd.merge(df,df2.loc[:,['Ensembl gene ID','Approved symbol']],how='left',on = 'Ensembl gene ID')

tcga=tcga.set_index(['Approved symbol'])

f=open("/media/liyaru/LYR/colon/27TCGA/TCGA/meta_tcga.txt")
T=[]
N=[]
lines=f.readlines()
for l in lines:
    l=l.split("\t")
    if l[1]=='Tumor':
        T.append(l[0])
    if l[1]=="Normal":
        N.append(l[0])
        
tcga["mean_N"]=tcga[N].mean(axis=1)
tcga["mean_T"]=tcga[T].mean(axis=1)
print(ttest_rel(tcga["mean_N"],tcga["mean_T"]))

down_SDOC_gene3=[]
for e in down_SDOC_gene2:
    if e in tcga.index:
        down_SDOC_gene3.append(e)
              
diff_SDOC_tcga=tcga.loc[down_SDOC_gene3,:]

diff_SDOC_tcga.to_csv(path_or_buf="./analysis_data/down_SDOC_TCGA_exp",sep="\t")
#diff_SDOC_tcga.to_csv(path_or_buf="./analysis_data/up_SDOC_TCGA_exp",sep="\t")

print(wilcoxon(diff_SDOC_tcga["mean_N"],diff_SDOC_tcga["mean_T"],alternative="less"))
print(wilcoxon(diff_SDOC_tcga["mean_N"],diff_SDOC_tcga["mean_T"],alternative="greater"))
print(ttest_rel(diff_SDOC_tcga["mean_N"],diff_SDOC_tcga["mean_T"]))
plt.figure()
f=plt.boxplot([diff_SDOC_tcga["mean_N"],diff_SDOC_tcga["mean_T"]],labels=["N","T"],widths=0.1,positions=[0,0.15])
