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
    outp="%s.windows%sk" % (args.o,args.s)
    out1=open("%s.site.txt"%(args.o) ,"w"  )
    out2=open("%s.txt"% (outp),"w" )
    bins=[]
    windows={}

    for line in open(args.v):
        if line.startswith("##"):
            continue
        elif line.startswith("#CHR"):
            name=line.strip().split()
            out1.write("Chr\tPOS\tRef\tAlt\ted2\ted4\n")
            out2.write("Chr\tPOS\tEnd\tNumber\tmeanED2\tmeanED4\n")
        else:
            info = line.strip().split()
            if "chr" not in info[0]: 
                continue
            ed2 = calED(info[9],info[10])
            if ed2 == "NN":
                continue
            ed4 = ed2**2
            out1.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (info[0],info[1],info[3],info[4],ed2,ed4))
            binID = int(info[1])/size 
            name = "%s\t%s" %(info[0],binID)
            if name not in windows:
                bins.append(name)
                windows[name] = {}
                windows[name]["ed2"] = []
                windows[name]["ed4"] = []
            windows[name]["ed2"].append(ed2)
            windows[name]["ed4"].append(ed4)

    out1.close()

    keep=0
    for name in bins: 
        chrID,binID = name.split("\t")
        if len(windows[name]["ed2"]) < args.min:
            continue
        keep += 1
        meanED2 = sum(windows[name]["ed2"])/len(windows[name]["ed2"])
        meanED4 = sum(windows[name]["ed4"])/len(windows[name]["ed4"])
        start = int(binID)*size+1 
        end   = (int(binID)+1)*size
        out2.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrID,start,end,len(windows[name]["ed4"]),meanED2,meanED4))
    out2.close()
    print "# Notice: At windows min site = %s, only %s out of %s were kept" % (args.min,keep,len(bins))

    cut1=draw("%s.site.txt"%(args.o),"%s.txt"% (outp),"ed2","meanED2","%s.ed2" % (outp))
    cut2=draw("%s.site.txt"%(args.o),"%s.txt"% (outp),"ed4","meanED4","%s.ed4" % (outp))
    with open("%s.cutline.txt" %(outp),"w") as outc:
        outc.write("#Cut-line\ned2\t%s\ned4\t%s\n" %(cut1,cut2))


def draw(pointfile,linefile,posC,lineC,outfile):
    ed=pd.read_csv(pointfile,sep='\t')
    windows= pd.read_csv(linefile,sep='\t')
    xmax=ed['POS'].max()
    ymax=ed[posC].max()

    cut_index =  int(0.05*len(windows[lineC]))
    cuth     =  sorted(windows[lineC],reverse=True)[cut_index]
    windows[windows[lineC] > cuth].to_csv( "%s.candidate.csv" % (outfile),sep='\t',index=False) 

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
        #axes[i].get_xaxis().get_major_formatter().set_useOffset(False)
    plt.savefig(fname="%s.png" % (outfile),dpi=300)
    plt.show()
    return cuth

def calED(s1,s2):
    a1,a2 = getAllen(s1.split(":")[1])
    b1,b2 = getAllen(s2.split(":")[1])
    total =  a1 + a2 + b1 + b2
    if  total < 10 or a1+a2 < 5 or b1+b2 < 5 :
         return "NN"
    ed = (float(a1)/(a1+a2)-float(a2)/(a1+a2))**2 + (float(b1)/(b1+b2)-float(b2)/(b1+b2))**2
    return ed 

def millions(x,pos):
    return '{:1.1f}M'.format(x*1e-6)

def getAllen(s):
    d1,d2 = s.split(",")
    return int(d1),int(d2)

def get_args():
    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton:  Performing Bulk-sequencing Analysis using Euclidean-distance.
    Writer:    Meng Minghui < mengminghui@whdiggers.com >
'''
    )
    parser.add_argument('-v',metavar='name',help='vcf file',required=True,type=str)
    parser.add_argument('-s',metavar='size',help='size of windows, kb',type=int,default=500)
    parser.add_argument('-o',metavar='out',help='pre-name of out files',type=str,default="BSA")
    parser.add_argument('-min',metavar='min',help='min site of a window',type=int,default=20)
    args = parser.parse_args()
    return args 

if __name__ == '__main__':
    main()

