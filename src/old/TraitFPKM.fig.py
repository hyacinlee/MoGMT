#!/usr/bin/env python3 
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.graphics.gofplots import qqline
import sys 


def main():
    if len(sys.argv) < 3:
        print ("%s  < fpkm file > < Phe file > < Gene  and Trait list >" % (sys.argv[0]))
        exit(1)
    
    ### reads Phe and FPKM data 
    T,Traits,Tsample = readData(sys.argv[2],"h")
    G,Gsample,Genes  = readData(sys.argv[1],"v")

    print ("## Total Traits: %s \n\n## Total Genes: %s\n\n## Run Association ... ..." % (len(Traits),len(Genes)))


    for l in open(sys.argv[3],"r"):
        (gname,tname) = l.strip().split()
        PlotCor(G[gname],T[tname],gname,tname)
    

def PlotCor(g,t,gname,tname):
    x=[]
    y=[]
    for sam in g.keys():
        if sam not in t or t[sam] == "NA":
            continue 
        else:
            x.append(float(g[sam]))
            y.append(float(t[sam]))

    #X = sm.add_constant(x)
    #model = sm.OLS(y,X).fit()
    #predictions = model.predict(X)  

    plt.style.use('seaborn')
    fig,ax = plt.subplots()
    #ax = plt.subplot(111)
    ax.scatter(x, y)
    #ax.plot(x, predictions, 'r-', linewidth=2)
    ax.set_xlabel(gname, fontsize=12)
    ax.set_ylabel(tname, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.grid(False)
    qqline(ax, "r", x, y)
    plt.savefig('%s_vs_%s.png' %(tname,gname))


def readData(infile,order):

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
                    P=accumDict(P,k1,k2,v) 
                if order == "v":  ##  V-line  as key
                    P=accumDict(P,k2,k1,v)

    return P,Hlist,Vlist


def accumDict(di,k1,k2,v):
    if k1 not in di:
        di[k1] = {}
    else:
        di[k1][k2] = v 

    return di



main()


