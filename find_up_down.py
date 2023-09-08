# -*- coding:utf-8 -*-
# @Author:PengChen
# @Date:2023/4/5
# @File:find_up_down.py
# @Email:pengchen2001@zju.edu.cn

import pandas as pd
import os
import argparse
import time

parser = argparse.ArgumentParser()

parser.add_argument("--debug", action="store_true", help="enable debug mode", default=False)
parser.add_argument("-i","--input", help="input file 2 name")
parser.add_argument("-o","--output", help="output file name")
parser.add_argument("-l","--log", help="log file name")
parser.add_argument("-s","--sep",type=int, help="how many sep print log",default=100000)

args = parser.parse_args()

if args.debug:
    print("Debug mode enabled")
    input_file = '../test2'
    output_file = 'output.txt'
    log_file='log'
    n_sep=10
else:
    input_file = args.input
    output_file = args.output
    log_file=args.log
    n_sep=int(args.sep)

def find_closest(num, data):
    # 使用二分查找算法查找最接近num的位置
    left, right = 0, len(data) - 1
    while left <= right:
        mid = (left + right) // 2
        if data[mid] <= num:
            left = mid + 1
        else:
            right = mid - 1
    
    # 返回最接近num的位置
    return right

def find_row_range(num, first_column, second_column):
    if len(first_column)<1:
        return ([None,None])
    
    # 使用二分查找来找到最接近num的数字
    right=find_closest(num, first_column)
    if right<0:
        return ([None,0])
      
    # 检查num是否小于第二列
    if num <= second_column[right]:
        return ([right,right])
    else:
        if right==len(first_column)-1:
            return ([right,None])
        return ([right,right+1])

    
# find_row_range(16,[1,10,25],[5,20,30])
# find_row_range(-1,[1,10,25],[5,20,30])
# find_row_range(40,[1,10,25],[5,20,30])
# find_row_range(25,[1,10,25],[5,20,30])
# find_row_range(22,[1,10,25],[5,20,30])


#df = df.sort_values(by=['seq','start'])
#print('sort done')

if(os.path.exists(output_file)):
    os.remove(output_file)
if(os.path.exists(log_file)):
    os.remove(log_file)

all_res=[]
cnt=0
query_seq=""
prefix=""

with open(input_file, 'r') as f:
    for line in f:
        # 去除字符串开头和结尾的空白字符，并按制表符分割为列表
        row = line.strip().split('\t')
        #spacer target middle position
        mid_pos=(int(row[2])+int(row[3]))/2
        
        #get the seq   
        if query_seq!=row[1]:
            #open the file, 取了几个前缀就是几。
            if row[1][:4]!=prefix:
                prefix=row[1][:4]
                print('doing '+prefix)
                with open(log_file, 'a') as log:
                    log.write('doing '+prefix+'\n')
                if os.path.exists(f"each_seq/{prefix}.txt"):
                    df = pd.read_csv(f"each_seq/{prefix}.txt", sep='\t', usecols=[0,1,2,3], header=None, names=['id', 'seq', 'start', 'end'])

            query_seq=row[1]
            seq_df=df.query("seq == @query_seq")
            start=list(seq_df['start'])
            end=list(seq_df['end'])
            
        #find genes
        
        i=find_row_range(mid_pos,start,end)
        
        res=[row[0],"NA","NA"]
        
        if not i[0] is None:
            res[1]=' '.join(seq_df.iloc[i[0]].astype(str).tolist())
        
        if not i[1] is None:
            res[2]=' '.join(seq_df.iloc[i[1]].astype(str).tolist())
        
        all_res.append('\t'.join(res))    
        
        cnt+=1
        if(cnt%n_sep==0):
            with open(log_file, 'a') as log:
                log.write(time.strftime("%m-%d-%H:%M:%S", time.localtime())+'===========\n')
                log.write(f'{cnt} sequences done\n')
            # 将字符串写入新文件
            with open(output_file, 'a') as f:
                for i in all_res:
                    f.write(str(i)+'\n')
            all_res=[]
            
with open(output_file, 'a') as f:
    for i in all_res:
        f.write(str(i)+'\n')

#https://blog.csdn.net/lamusique/article/details/112430216?ydreferer=aHR0cHM6Ly93d3cuYmluZy5jb20v
#优化内存

#!/bin/bash
#SBATCH --job-name=find_up
#SBATCH --output=/data/home/jianglab/share/bac2pc/%x_%A.out
#SBATCH --error=/data/home/jianglab/share/bac2pc/%x_%A.err
#SBATCH --time=14-00:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

# echo start: `date +'%Y-%m-%d %T'`
# start=`date +%s`
# 
# ./find_up_down.py -i 3bac_total_bac_in_gene -o spacer_up_down
# 
# echo end: `date +'%Y-%m-%d %T'`
# end=`date +%s`
# echo TIME:`expr $end - $start`s

# python3 find_up_down.py -i test2 -o spacer_up_down
