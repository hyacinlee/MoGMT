#!/usr/bin/env python3 
import statsmodels.api as sm
import sys 
import matplotlib.pyplot as plt
from statsmodels.graphics.gofplots import qqline
from scipy.stats import pearsonr
import Common
import argparse

def get_args(args):
    parser = argparse.ArgumentParser(description="Association analysis between gene expression matrix and phenotype in population",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-e", "--express",required=True, type=str,help="Path of the input gene expression matrix file")
    parser.add_argument("-p", "--phe",required=True, type=str,help="Path of the input phenotype file, *phe")
    parser.add_argument("-t", "--trait",required=True, type=str,help="Trait Name")
    parser.add_argument("-l", "--loci",help="Path of the gene loci file in bed format with 5 Colums")
    parser.add_argument("-o", "--out",type=str,help="Out frefix",default="output")
    parser.add_argument("--figGenes",type=str, help="a list file, Out Fpkm - Phe correlation plot with figGenes list file")
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args



def main(args=None):
    args=get_args(args)

    out=args.out    
    ### reads Phe 
    #P,Traits,Tsample = readData(args.phe,"h")
    Phe=Common.read_file(args.phe,mode="dict",vals=[2],keys=[0])
    ### reads FPKM data 
    G,Gsample,Genes  = readData(args.express,"v")


    print ("## Total Genes: %s\n\n## Run Association ... ..." % (len(Genes)))

    cut = 1.0/(len(Genes))
    with open(f"./Associate.cutline.txt","w") as cf:
        cf.write(str(cut))

    if args.figGenes:
        fGenes =  Common.read_file(args.figGenes,mode="list",vals=[0])
        for gene in fGenes:
            if gene not in Genes:
                continue
            single_lineRegress(Phe,G[gene],outplot=True,gname=gene,tname=out)
    else:
        run_lineRegress(args.trait,Phe,G,Genes,cut,"TransAsso.%s" %(out),args.loci)

  

def run_lineRegress(trait,Phe,Ginfo,Genes,cut,outfile,loci):
    
    Pos={}
    ### read Pos 
    if loci:
        of2=open(outfile+".rdata","w")
        of2.write("SNP\tChromosome\tPosition\tPvalue\n")

        for l in open(loci,"r"):
            ls=l.strip().split()
            mid = int((int(ls[1])+int(ls[2]))/2)
            Pos[ls[4]] = [ls[0],mid]


    of1=open(outfile+".txt","w")
    of1.write("Gene\tTrait\tlm_r2\tlm_b\tlm_p\tpearson_r\tpearson_p\n")
        
    of3=open(outfile+".signal.txt","w")
    of3.write("Gene\tTrait\tlm_r2\tlm_b\tlm_p\tpearson_r\tpearson_p\n")
    
    for gene in Genes:
        model,pearson = single_lineRegress(Phe,Ginfo[gene])
        r_squared = model.rsquared
        b = model.params[1]
        p_value = model.pvalues[1] if not str(model.pvalues[1]) == "nan" else  1.0

        #print(p_value)
        #print(cut)
        of1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,trait,r_squared,b,p_value,pearson[0],pearson[1]))

        if float(p_value) <= cut:
            of3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,trait,r_squared,b,p_value,pearson[0],pearson[1]))

        if loci:
            of2.write("%s\t%s\t%s\t%s\n" % (gene,Pos[gene][0],Pos[gene][1],p_value))


    of1.close()
    of3.close()
    if loci:
        of2.close()


def single_lineRegress(Phe,Express,outplot=False,gname=None,tname=None):
    x=[]
    y=[]
    for sam in Phe.keys():
        if Phe[sam] == "NA" or sam not in Express:
            continue
        else:
            x.append(float(Express[sam]))
            y.append(float(Phe[sam]))
              
    #print(x)
    xx = sm.add_constant(x)
    model = sm.OLS(y,xx).fit()

    pearson = pearsonr(x, y)

    if not outplot:
        return model,pearson

    

    r_squared = model.rsquared 
    p_value = model.pvalues[1] 

    #print(len(x),len(y))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x, y)
    ax.set_xlabel(gname, fontsize=12)
    ax.set_ylabel(tname, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)

    # 画拟合线
    ax.plot(x, model.fittedvalues, color='red')

    # 在图上标注 R² 和 P
    ax.text(
        0.05, 0.95,
        f"R² = {r_squared:.3f}\nP = {p_value:.3e}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)
    )
    plt.savefig('TransAsso.%s_vs_%s.png' %(tname,gname),dpi=300, bbox_inches='tight')



def single_cor_plot(Phe,Express,gname,tname):
    x=[]
    y=[]
    for sam in Phe.keys():
        if sam not in Express or Phe[sam] == "NA":
            continue 
        else:
            x.append(float(Express[sam]))
            y.append(float(Phe[sam]))

    x = sm.add_constant(x)
    model = sm.OLS(y,x).fit()
    r_squared = model.rsquared 
    p_value = model.pvalues[1] 


    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x, y)
    ax.set_xlabel(gname, fontsize=12)
    ax.set_ylabel(tname, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)

# 画拟合线
    ax.plot(x, model.fittedvalues, color='red')

# 在图上标注 R² 和 P
    ax.text(
        0.05, 0.95,
        f"R² = {r_squared:.3f}\nP = {p_value:.3e}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)
    )
    plt.savefig('TransAsso.%s_vs_%s.png' %(tname,gname),dpi=300, bbox_inches='tight')








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


if __name__ == "__main__":
    main()
