#!/usr/bin/env python3 
import sys 
import Common
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from qmplot import manhattanplot



def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-i","--input",required=True, help="Gwas or Other site-value Result")
    parser.add_argument("-c","--chrom",help="Chromosom name column",type=str,default="CHROM") 
    parser.add_argument("-p","--pos",help="Maker pos column",type=str,default="POS")
    parser.add_argument("-v","--value",help="Value column",type=str,default="P")
    parser.add_argument("-n","--name",help="SNP name column",type=str,default="ID")
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out")
    parser.add_argument("-s","--sign",help="Significance cut line: B = 1/n || F = 0.05/n || top1 = top 0.01 || top5 = top 0.05  || a self-defined float",default="B")
    parser.add_argument("-l","--log",help="change to log10 value",action="store_true",default=True)
    parser.add_argument("--sub",help="Draw only one chr",type=str,default=None)
    parser.add_argument("--hightlight",help="Hightlight maker ID in files, ids must be same with input",type=str,default=None)


    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args
 

def main(args=None):
 
    args=get_args(args)
    ManhantanMain(input=args.input,chrom=args.chrom,pos=args.pos,value=args.value,name=args.name,out=args.out,sign=args.sign,log=args.log,sub=args.sub,hightlight=args.hightlight)


def ManhantanMain(input,chrom="CHROM",pos="POS",value="P",name="ID",out="out",sign="B",log=True,sub=None,hightlight=None):
    
    name_col= {chrom:"#CHROM",pos:"POS",value:"P",name:"ID"}
    df = pd.read_table(input, sep="\t")
    df = df.rename(columns=name_col)
    df = df.dropna(how="any", axis=0)

    #print(sign)
    cut = Common.judge_significant(sign,df)
 
    #filtered_df = df[df['P'] < cut]
    #filtered_df.to_csv("%s.sign.xls" %(o),sep="\t",index=False)

    yname="-log10(P)"
    cut = -np.log10(cut)
    if not log:
        df["P"]=1.0/10**df["P"]
        yname=value
        cut = 1/10**float(cut)

    manhatan_fig(df,cut,yname,out=out,sub=sub,hightlight=hightlight)
    print(f"# Finish ouptut manhantun plot in {out}.png")
    #ManhantunMain(input,chrom=args.chrom,pos=args.pos,value=args.value,out=args.out,sign=args.sub,log=args.log,sub=args.sub,hightlight=args.hightlight)


def manhatan_fig(df,cut,ylab="-log10(P)",out="out",sub=None,hightlight=None):
    #print df
    ax=""
    out="%s.png" %(out)
    if not sub:
        ax=manhattanplot(data=df,xlabel="Chromosome",ylabel=ylab,suggestiveline=None,genomewideline=None,) 
    else:
        #out="%s.%s.pdf" %(args.o,args.sub)
        ax=manhattanplot(data=df,xlabel=sub,CHR=sub,ylabel=ylab,suggestiveline=None,genomewideline=None)
    
    ax.axhline(y=cut, linewidth = 2,linestyle="--",color="red")
    plt.savefig(out,orientation="landscape",bbox_inches='tight',dpi=300)
    plt.close()



def significanceDF(df,cut,out):
    sites={}
    GwasP={}
    sitesl=[]

    signDf=df[df['P'] < cut ]
    # out signDf  
    logger.info("%s significance SNP were founded, wirting them into %s.sign.sites"%(len(signDf),out))
    signDf.to_csv("%s.sign.sites" %(out),index=False,sep="\t")
 
    for index,row in signDf.iterrows():
        c=str(row['#CHROM'])
        p=int(row['POS'])
        i=str(row['ID'])
        v=float(row['P'])
        sites=accumulateDict(sites,i,c,p) 
        sitesl.append([i,c,p])
        GwasP[i] = v

    return (sites,sitesl,GwasP)



if __name__ == '__main__':

    main()
