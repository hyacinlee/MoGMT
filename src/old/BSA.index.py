#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():
    args=get_args()

    size=int(args.s)*1000
    out1=open("%s.site.txt"%(args.o) ,"w"  )
    out2=open("%s.%sk.windows.txt"% (args.o,args.s),"w" )
    bins=[]
    sample=[]
    windows={}

    for line in open(args.v):
        if line.startswith("##"):
            continue
        elif line.startswith("#CHR"):
            info=line.strip().split()
            sample = (info[9],info[10],info[11],info[12])
            out1.write("Chr\tPOS\tParent1\tParent2\tBulk1-depth\tBulk2-depth\tBulk1-index\tBulk2-index\tdeltaIndex\n")
            out2.write("Chr\tPOS\tEnd\tNumber\tWindeltaIndex\n")
        else:
            info = line.strip().split()
            gts = (info[9],info[10],info[11],info[12])
            ginf = getGT(sample,gts)
            #print ginf
            excp=["0/1","./.","1/0"]
            if ginf[args.p1]["GT"] == ginf[args.p2]["GT"] or ginf[args.p2]["GT"] in excp or ginf[args.p1]["GT"] in excp:
                continue
            if ginf[args.b1]["GT"] == "./." or ginf[args.b2]["GT"] == "./.":
                continue
            if ginf[args.b1]["D0"] + ginf[args.b1]["D1"] == 0 or ginf[args.b2]["D0"] + ginf[args.b2]["D1"] == 0:
                continue

            index1 = calIndex(ginf[args.p1]["GT"],ginf[args.b1])
            index2 = calIndex(ginf[args.p1]["GT"],ginf[args.b2])
            if index1 == 999 or index2==999:
                 print ginf
                 continue
            delta=index1-index2
            
            gp1 = getAllenSingle(info[3],info[4],ginf[args.p1]["GT"])
            gp2 = getAllenSingle(info[3],info[4],ginf[args.p2]["GT"])

            bg1 = "%s:%s|%s:%s" % (info[3],ginf[args.b1]["D0"],info[4],ginf[args.b1]["D1"])
            bg2 = "%s:%s|%s:%s" % (info[3],ginf[args.b2]["D0"],info[4],ginf[args.b2]["D1"])
            out1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (info[0],info[1],gp1,gp2,bg1,bg2,index1,index2,delta))
            binID = int(info[1])/size 
            name = "%s\t%s" %(info[0],binID)
            if name not in windows:
                bins.append(name)
                windows[name] = {}
                windows[name]["delta"] = []
            windows[name]["delta"].append(delta)

    out1.close()

    keep=0
    for name in bins: 
        chrID,binID = name.split("\t")
        if len(windows[name]["delta"]) < args.min:
            continue
        keep += 1
        meandelta = sum(windows[name]["delta"])/len(windows[name]["delta"])
        start = int(binID)*size+1 
        end   = (int(binID)+1)*size
        out2.write("%s\t%s\t%s\t%s\t%s\n" % (chrID,start,end,len(windows[name]["delta"]),meandelta))
    out2.close()
    print "# Notice: At windows min site = %s, only %s out of %s were kept" % (args.min,keep,len(bins))

    draw("%s.site.txt"%(args.o),"%s.%sk.windows.txt"% (args.o,args.s),"deltaIndex","WindeltaIndex","%s.windos%sk.delta" % (args.o,args.s),args.c)
    #with open("cutline.txt","w") as outc:
 #       outc.write("#Cut-line\ned2\t%s\n" %(cut1))
#

def getGT(name,dat):
    inf={}
    for i in range(len(name)):
        infs=dat[i].split(":")
        d =  infs[1].split(",")
        inf[name[i]]={}
        inf[name[i]]["GT"]=infs[0]
        inf[name[i]]["D0"]=int(d[0])
        inf[name[i]]["D1"]=int(d[1])
    return inf

def getAllen(ref,alt,gt):
    if gt =="0/0":
        return ref+ref  
    elif gt == "0/1":
        return ref+alt
    elif gt == "1/1":
        return alt+alt
    else:
        return "NN"

def getAllenSingle(ref,alt,gt):
    if gt =="0/0":
	return ref
    elif gt == "1/1":
        return alt
    else:
        return "N"


def draw(pointfile,linefile,posC,lineC,outfile,cut):
    ed=pd.read_csv(pointfile,sep='\t')
    windows= pd.read_csv(linefile,sep='\t')
    xmax=ed['POS'].max()
    ymax=ed[posC].max()

    #cut_index =  int(0.05*len(windows[lineC]))
    #cuth     =  sorted(windows[lineC],reverse=True)[cut_index]
    cuth = cut
    cutl = 0.0-cut
    windows[ (windows[lineC] > cuth) | (windows[lineC] < cutl)].to_csv( "%s.candidate.csv" % (outfile),sep='\t',index=False) 

    bins=list(windows['Chr'].unique())
    fig,axes = plt.subplots(len(bins),1,sharex=True,figsize=(2*len(bins),len(bins)))
    for i in range(len(bins)):
        name = bins[i]
        axes[i].text(xmax*1.05,ymax/2,name,fontsize=15,verticalalignment="top",horizontalalignment="left")
        axes[i].scatter(ed[ed['Chr']==name]['POS'],ed[ed['Chr']==name][posC],s=0.5,c='lightgray')
        axes[i].plot(windows[windows['Chr']==name]['POS'],windows[windows['Chr']==name][lineC],c="firebrick")
        axes[i].spines['top'].set_color("white")
        axes[i].spines['right'].set_color("white")
        axes[i].xaxis.set_major_formatter(plt.FuncFormatter(millions))
        axes[i].axhline(y=cuth,ls="--",c="navy",lw=1)
        axes[i].axhline(y=cutl,ls="--",c="navy",lw=1)
        #axes[i].get_xaxis().get_major_formatter().set_useOffset(False)
    plt.savefig(fname="%s.png" % (outfile),dpi=300)
    plt.show()
    return cuth


def calIndex(refGT,bulk):
    index=999
    if refGT == "0/0":
        index = bulk["D0"]/float((bulk["D0"]+bulk["D1"]))
    elif refGT == "1/1":    
        index = bulk["D1"]/float((bulk["D0"]+bulk["D1"]))
    return float(index)


def millions(x,pos):
    return '{:1.1f}M'.format(x*1e-6)

def get_args():
    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton:  Performing Bulk-sequencing Analysis using SNP-index.
    Writer:    Meng Minghui < mengminghui@whdiggers.com >
'''
    )
    parser.add_argument('-v',metavar='name',help='vcf file',required=True,type=str)
    parser.add_argument('-s',metavar='size',help='size of windows, kb',type=int,default=500)
    parser.add_argument('-p1',metavar='p1',help='name of Parent 1,as ref ',type=str)
    parser.add_argument('-p2',metavar='p2',help='name of Parent 2',type=str)
    parser.add_argument('-b1',metavar='b1',help='name of bulk1',type=str)
    parser.add_argument('-b2',metavar='b2',help='name of bulk2',type=str)
    parser.add_argument('-o',metavar='out',help='pre-name of out files',type=str,default="BSA")
    parser.add_argument('-min',metavar='min',help='min site of a window',type=int,default=20)
    parser.add_argument('-c',metavar='c',help='cut off',type=float,default=0.5)
    args = parser.parse_args()
    return args 

if __name__ == '__main__':
    main()
