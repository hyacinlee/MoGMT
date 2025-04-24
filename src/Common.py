#!/usr/bin/env python3 
import re
import os
import sys 
import argparse
import subprocess
import logging
import numpy as np
import pandas as pd
from scipy import stats


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawTextHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        return ', '.join(action.option_strings)

    def _get_help_string(self, action):
        help_str = super()._get_help_string(action)
        
        if action.required:
            help_str += ' (required)'

        #if action.default:
        #    help_str = help_str.replace("(default: None)", "")
        
        return help_str



# File - List & Dict read/write  Function 

def read_file(infile,mode="list",vals=[],keys=[],header=False,sep="\t"):
    result = [] if mode == "list" else {}
    head_info=[]
    with open(infile,"r") as inf:
        for i, line in enumerate(inf):
            if line.startswith("#"):    # start with #
                continue

            datas = line.strip()
            if not sep == None:   
                ss=line.strip().split(sep)
                datas = ss if len(vals) == 0 else [ ss[i] for i in vals ]

            if header and i==0:
                head_info=datas
                continue

            ss = line.strip().split(sep)
            #datas = ss if len(vals) == 0 else [ ss[i] for i in vals ]

            if mode == "list":
                result.append(datas)
            elif mode == "dict":
                key = "###".join([ss[i] for i in keys])
                result[key] = datas

    if header==False:
        return result
    else:
        return (result,head_info)


def read_file_accumulateDict(infile,vals=[0],key1=1,key2="no",sep="\t"):

    result={}
    with open(infile,"r") as inf:
        for line in inf:
            if line.startswith("#"):    # start with #
                continue

            datas = line.strip().split(sep) 
            vv=""
            if len(vals)==1:  # only one 
                vv=datas[vals[0]]
            else:
                vv=[datas[i] for i in vals] 
  
            if not key2 =="no":
                key2 = datas[key2]

            result=accumulateDict(result,vv,datas[key1],key2)
            
        return result












def write_list_to_file(mylist,output,vals=None,header=False,sep="\t"):
    # write a list to file
    # vals = None : a simple list 
    # vals != None: a two dimom list

    with open(output, 'w') as f:

        if vals == None:
            if header:
                f.write(header+"\n")
            for s in mylist:
                f.write(s+"\n")
        else:
            if header:
                datas = header if len(vals) == 0 else [ header[i] for i in vals ]
                f.write(sep.join(datas)+"\n")
            for row in mylist:
                row = list_type(row,"str")
                datas = row if len(vals) == 0 else [ row[i] for i in vals ]
                f.write(sep.join(datas)+"\n")


def sort_dict_by_list(mydict,sort_list):
    # sort and filter a dict by a sorted list
    # return the sorted dict value and list 
    vs=[]
    ks=[]
    for s in sort_list:
        if s in mydict:
            ks.append(s)
            vs.append(mydict[s])
    return vs,ks


def read_matrix_data(infile,order):
    P={}
    Hlist=[]
    Vlist=[]
    for l in open(infile,"r"):
        if "t_name" in l:
            continue

        ls = l.strip().split()
        name=ls.pop(0)

        if len(Hlist) ==0:
            Hlist=ls
        else:
            Vlist.append(name)
            for i in range(len(Hlist)):
                k1 = Hlist[i]
                k2 = name
                v = ls[i]
                if order == "h":  ##  H-line  as key
                    P=accumulateDict(P,v,k1,k2) 
                if order == "v":  ##  V-line  as key
                    P=accumulateDict(P,v,k2,k1)
    return P,Hlist,Vlist


def accumulateList(mylist,val):
    if val not in mylist:
        mylist.append(val)
    return mylist


def accumulateDict(myDict,val,keyA,keyB="no"):
    '''
        accumulate a 1~2d Dictionary
    '''
    if keyB == "no":
        if keyA in myDict:
            if type(val) == int :
                myDict[keyA] += val
            else:
                myDict[keyA].append(val)         
        else:
            if type(val) == int :
                myDict[keyA] = val
            else:
                myDict[keyA] = [val]

    else:
        if keyA in myDict:
            if keyB in myDict[keyA]:
                myDict[keyA][keyB] += val
            else:
                myDict[keyA].update({keyB: val})   
        else:
            myDict.update({keyA:{keyB: val}})

    return myDict


def fillDictValue(mydict,keylist,val):
    for k in keylist:
        if k not in mydict:
            mydict[k] = val
    return mydict


def list_type(mylist,types):
    nl=[]
    for s in mylist:
        if types == "str":
            nl.append(str(s))
        elif types == "float":
            nl.append(float(s))
    return nl


# Syetem Function 
def run_command(cmd, check=True):
    """Helper function to run shell commands"""
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=check)


def check_path_exists(path):
    if not os.path.exists(path):
        print(f"Error：Paht of software or file not exist: {path}")
        sys.exit(1) 
    else:
        return True


def configure_logging():
    #logging.basicConfig(filename="EGCV.log",level=logging.INFO,format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",)
    return logging.getLogger("EVCG")



## Data processing Function 
def judge_significant(sign,df):
    cut = sign 

    if cut == "B":
        cut=float(1)/len(df)
        print(f"# The significance threshold is {cut} ")
    elif cut == "F":
        cut=float(0.05)/len(df)
        print(f"# The significance threshold is {cut} ")
    else:
        cut = float(cut)

    return cut


def search_cloesd_region(gene,locis):
    # gene is a list with [chrom,start,end,strand] 
    # loci is a lsit with [[chrom,start,end,id],...]

    #print(gene)
    #print(locis[0:10])

    chrom_g, start_g, end_g, strand_g = gene[0:4]
    overlapping_regions = []
    
    for locus in locis:
        chrom_l, start_l, end_l, id_l = locus[0:4]
        
        if chrom_l != chrom_g:
            continue
        
        if (int(start_l) <= int(end_g)) and (int(end_l) >= int(start_g)):
            overlapping_regions.append(locus)

    return overlapping_regions


def rename_dataframe(df,chrom,pos,name,value):
    name_col= {chrom:"#CHROM",pos:"POS",value:"P",name:"ID"}
    df = df.rename(columns=name_col)
    return df



if __name__ == '__main__':
    logger = configure_logging()
    main()
