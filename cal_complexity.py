# -*- coding:utf-8 -*-
# @Author:PengChen
# @Date:2022/11/23
# @File:cal_complexity.py
# @Email:pengchen2001@zju.edu.cn

import Bio.SeqUtils.lcc
from Bio import SeqIO
import argparse,os,sys,math
from collections import Counter
import pandas as pd

parser = argparse.ArgumentParser(description='This script is used to calculate sequence complexity')
parser.add_argument('-s','--string',help='a DNA/RNA sequence',required=False)
parser.add_argument('-f','--file',help='Please input a fasta file',required=False)
parser.add_argument('-o','--out_dir',help='Please input out_put directory path',default = os.getcwd(),required=False)

args = parser.parse_args()
# seq="ATCGTGGCAGT"
# wsize=3
# #similar to shannon
# #Bio.SeqUtils.lcc.lcc_simp("ATC")
# complexity=Bio.SeqUtils.lcc.lcc_mult(seq, wsize)
# import matplotlib.pyplot as plt
# site=list(range(1,len(seq)+2-wsize))
# plt.plot(site, complexity, 'b*--',label='Complexity')
# plt.legend()
def count_kmers(sequence, k_size):
    data = {}
    size = len(sequence)
    for i in range(size - k_size + 1):
        kmer = sequence[i: i + k_size]
        try:
            data[kmer] += 1
        except KeyError:
            data[kmer] = 1
    return data

def fenmu(i,lens):
    a=4**i
    b=lens+1-i
    if(a>b):
        a=b
    return a

def cal_sc(sequence,max_k=None):
    sc=1
    L=len(sequence)
    if max_k is None:
        max_k=7
    max_k=min(max_k,L)
    for k in range(1,max_k):
        sc=sc*len(count_kmers(sequence,k))/fenmu(k,L)
    return round(sc,5)

def summary_seq(name,sequence,max_k=7):
    global counts
    counts["id"].append(name)
    counts["length"].append(len(sequence))
    tmp = Counter(sequence)
    for i in ["A", "T", "C", "G"]:
        counts[i].append(tmp[i])
    counts["shannon"].append(round(Bio.SeqUtils.lcc.lcc_simp(sequence),5))
    counts["complexity"].append(cal_sc(sequence,max_k))

counts= {"id":[],"shannon":[],"complexity":[],"length":[],"A":[],"T":[],"C":[],"G":[]}
if (args.string):
    summary_seq("your sequence",args.string)
    print(counts)
    exit()
if (args.file):
    file=args.file
    for record in SeqIO.parse(file, "fasta"):
        summary_seq(record.name,record.seq)
    res=pd.DataFrame(counts)
    out=args.out_dir+"/"+os.path.basename(file).split(".")[0]+"_complexity.csv"
    res.to_csv(out,index=False)
    print("done! see " + out)



#=====max_k选7就够了
if (False):
    file="test.fa"
    df=pd.DataFrame(dict.fromkeys(range(3,12),[]))
    for record in SeqIO.parse(file, "fasta"):
        tmp=[cal_sc(record.seq,i) for i in range(3,12)]
        df.loc[df.index.size]=tmp
    x=range(3,12)
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use('TkAgg')
    for i in range(0,5):
        plt.plot(x, df.loc[i], 'r', marker='*', markersize=10)
    plt.show(block=True)
    plt.close()

#=====随机取序列碎片
if(False):
    import random
    import numpy as np

    with open("random_seqs.fa", "w") as f:
        for i in range(10000):
            f.write(">random" + str(i) + "\n")
            seq = random.choices(["A", "T", "C", "G"], k=int(round(np.random.normal(loc=33,scale=2),0)))
            seq = "".join(seq)
            f.write(seq + "\n")

#======真实bacteria genome repeats取样
if(False):
    import random
    import numpy as np
    import pandas as pd
    def rep_pick(r_name,r_seq):
        n=int(len(r_seq) / 10)
        strs=""
        for i in range(5):
            start = random.randint(0, len(r_seq) - 33)
            k = int(round(np.random.normal(loc=33, scale=2), 0))
            tmp_seq = r_seq[start:start + k]
            strs=strs+">" + r_name + "_" + str(start) + ":" + str(start + k - 1)+"\n"+tmp_seq+"\n"
        return strs

    #reps=pd.read_csv("/Users/asa/Documents/R/crisper/data/rep_10000.tsv",sep="\t")
    reps=pd.read_csv("/Users/asa/Documents/R/crisper/data/reps_20.tsv",sep="\t")
    with open("/Users/asa/Documents/R/crisper/data/bacteria_tandem_repeat20.fa","w")as f:
        for i in range(reps.shape[0]):
            f.write(rep_pick("Row"+str(i),reps.loc[i,"sequence"]))



