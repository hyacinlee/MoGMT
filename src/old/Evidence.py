#!/usr/bin/env python3 
import re
import sys 
import logging
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from qmplot import manhattanplot
from scipy import stats
from cyvcf2 import VCF

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass

class NoSubparsersHelpFormatter(argparse.HelpFormatter):
    def add_arguments(self, actions):
        actions = [a for a in actions if not isinstance(a,argparse._SubParsersAction)]
        super().add_arguments(actions)

def get_args():
    # parent_parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    # common setting 
    parent_parser.add_argument("-h", "--help", action="help", help="help")
    parent_parser.add_argument('-i',metavar='',help='Gwas or Other site-value Result',type=str)
    parent_parser.add_argument('-chr',metavar='',help='chromosom  name',type=str,default="Chr")
    parent_parser.add_argument('-pos',metavar='',help='maker pos name',type=str,default="Pos")
    parent_parser.add_argument('-value',metavar='',help='value name',type=str,default="Pvalue")
    parent_parser.add_argument('-log',help='change to log10 value',action="store_true",default=False)
    parent_parser.add_argument('-o',metavar='',help='out prefix',default="out")

    # main_parser
    parser = argparse.ArgumentParser(
        formatter_class = HelpFormatter,
        description = "RunCmd:\n\tFunciton:\n\tWriter: Meng Minghui <hyacinlee@163.com>")
    subparsers = parser.add_subparsers(dest="command")

    # Evidence parser
    parser_Evidence = subparsers.add_parser("Evidence",parents=[parent_parser],add_help=False,
        formatter_class=NoSubparsersHelpFormatter,
        help="Evidence GWAS result in Genes",
        description="Evidence GWAS result in Genes")
    parser_Evidence.add_argument('-phe',metavar='',help='Phenotype table of each sample in three cols: indv | indv | value',type=str)  
    #parser_Evidence.add_argument('-express',metavar='',help='express: fpkm of mRNA',type=str)
    parser_Evidence.add_argument('-vcf',metavar='',help='vcf file in bgzip format and has bcftools-csi index',type=str)
    parser_Evidence.add_argument('-annovar',metavar='',help='snp annovar file, *variant_function',type=str)
    parser_Evidence.add_argument('-annovar_exonic',metavar='',help='snp annovar exonic file, *exonic_variant_function',type=str)
    parser_Evidence.add_argument('-weight',metavar='',help='Other evidence with weight',type=str)

    # manhatan parser
    parser_Manhantun = subparsers.add_parser("Manhantun",parents=[parent_parser],add_help=False,
        formatter_class=NoSubparsersHelpFormatter,
        help="Plotting Gonome-wide or sub-Chr manhatan plot use any Dataset",
        description="Evidence GWAS result in Genes")
    parser_Manhantun.add_argument('-sub',metavar='',help='draw only one chr')
    parser_Manhantun.add_argument('-sign',metavar='',help='significance cut line: B = 1/n || F = 0.05/n || top1 = top 0.01 || top5 = top 0.05 or || a self-defined float',type=str,default="B")
    parser_Manhantun.add_argument('-hightlight',metavar='',help='hightlight maker ID, must be same with input',type=str)

    
    parser_Twas = subparsers.add_parser("TWAS", help="Preform standalone FPKM--Phe Association analysis without Genotypes",
        description="Preform standalone FPKM--Phe Association analysis without Genotypes",
        parents=[parent_parser],add_help=False,formatter_class=NoSubparsersHelpFormatter)
    parser_Twas.add_argument('-phe',metavar='',help='Phenotype table of each sample in three cols: indv | indv | value',type=str)
    parser_Twas.add_argument('-express',metavar='',help='express: fpkm of mRNA',type=str)

    parser_PHEcluster = subparsers.add_parser("PHEcluster", help="Cluster Phenotype data with PCA and valued in PC1",
        description="Cluster Phenotype data with PCA and valued in PC1",
        parents=[parent_parser],add_help=False,formatter_class=NoSubparsersHelpFormatter)
    parser_PHEcluster.add_argument('-phe',metavar='',help='Phenotype table of each sample in three cols: indv | indv | value',type=str)
    parser_PHEcluster.add_argument("-b", type=int, help="xxxx")


    args = parser.parse_args()

    if args.command == None:
        parser.print_help()
        print ("\n### Please choose a sub function !\n")
        exit (1)

    return args   


def configure_logging():
    #logging.basicConfig(filename="EGCV.log",level=logging.INFO,format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",)
    return logging.getLogger("EVCG")


def main(): 
    args=get_args()
    logger.info('Starting main Program')
    logger.info(args)
    if args.command == "Evidence":
        EvidenceMain(args)
    elif args.command == "Manhantun":
        ManhantunMain(args)
    else:
        print ("\n### Please choose a sub function!\nUse %s -h/--help for more help infomation\n" % sys.args[0])


def ManhantunMain(args):

    name_col= {args.chr:"#CHROM",args.pos:"POS",args.value:"P",args.i:"ID"}

    df = pd.read_table(args.i, sep="\t")
    df = df.rename(columns=name_col)
    df = df.dropna(how="any", axis=0)

    cut =judgeCut(args.sign,df)
 
    filtered_df = df[df['P'] < cut]
    filtered_df.to_csv("%s.sign.xls" %(args.o),sep="\t",index=False)

    yname="-log10(P)"
    cut = -np.log10(cut)
    if not args.log:
        df["P"]=1.0/10**df["P"]
        yname=args.value
        cut = 1/10**float(cut)

    manhatan_fig(df,cut,yname,args)


def EvidenceMain(args):
    out= args.o

    logger.info("Reading GWAS result from %s" %(args.i))
    df=pd.read_table(args.i, sep="\t")
    cut=float(1)/len(df)  
    logger.info("Gwas cutline is %s"%(cut))
    (sites,sitesl,GwasP) = significanceDF(df,cut,out)

    logger.info("Reading Candidata SNP Genotype info from %s" %(args.vcf))
    (cgt,cinfo)=readSNPs(args.vcf,sitesl,out)

    logger.info("Reading significance SNP function from annovar result %s and %s" %(args.annovar,args.annovar_exonic)) 
    (SNPcount,SNPloci,snp2gene) = readAnnovar(args.annovar,args.annovar_exonic,sites) 

    logger.info("Reading Phenotype data from %s" %(args.phe)) 
    Phe = readTableDict(args.phe,0,2)

    (snp,gene)=outBasicModel(out,sites,GwasP,cgt,cinfo,SNPcount,SNPloci,snp2gene,Phe) 
    
    if args.weight:  # with eight model
        outWeightModel(args.weight)

def outWeightModel(weightFile):
    evidens=readData(weightFile,"v")
    print (evidens)


        #logger.info("Reading the gene expression matrix from %s" %(args.express))       
        #G,Gsample,Genes  = readData(args.express,"v")
        #pass 
        #    for recode in snp2gene[snp]:         # student t-test SNP vs Express
        #        g=recode[0]
        #        if g == "NONE":
        #            logger.warning("snp %s in Gene %s has none express data, skipping" %(snp,g))
        #            continue
        #        (pp2,p2,r)=stTest(cgt[sID],G[g])         
        #if g not in G:
        #    logger.warning("%s has no Transcript Data, skipping" % (g))
        #else: 
        #    lineregress=expressPheAsso(G[g],Phe)
        #    p=lineregress.pvalues[1]
    #manhatan(sys.argv[1],"%s.png" %(sys.argv[1]))

def manhatan_fig(df,cut,yname,args):
    #print df
    ax=""
    out="%s.png" %(args.o)
    if not args.sub:
        ax=manhattanplot(data=df,xlabel="Chromosome",ylabel=yname,suggestiveline=None,genomewideline=None,) 
    else:
        #out="%s.%s.pdf" %(args.o,args.sub)
        ax=manhattanplot(data=df,xlabel=args.sub,CHR=args.sub,ylabel=yname,suggestiveline=None,genomewideline=None)
    
    ax.axhline(y=cut, linewidth = 2,linestyle="--",color="red")
    plt.savefig(out,orientation="landscape",bbox_inches='tight',dpi=300)


def outBasicModel(out,sites,GwasP,cgt,cinfo,SNPcount,SNPloci,snp2gene,Phe):
    snp=[]
    gene=[]
    logger.info("Outputing significance SNP function, GWAS-Pvalue")
    with open("%s.sign.BasicSite.xls" % (out),"w") as of1:
        of1.write("#ID\tRef\tAlt\tType\tGene\tP-GWAS\tP-CorPhe\n")
        for snp in snp2gene.keys():
            rs = snp.split("#")
            if rs[0] not in sites or int(rs[1]) not in sites[rs[0]]:
                logger.warning("Check SNP %s %s" %(rs[0],rs[1]))
                continue
            sID = sites[rs[0]][int(rs[1])]
            (pp1,p1,r)=stTest(cgt[sID],Phe)     # student t-test SNP vs Phe
            for recode in snp2gene[snp]:
                of1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(sID,cinfo[sID][0],cinfo[sID][1][0],recode[0],recode[1],GwasP[sID],p1))

    logger.info("Outputing Genes-SNP statistical table")
    with open("%s.sign.BasicGene.xls" % (out),"w") as of2:
        of2.write("#ID\tTotalSNP\tNoSys\tExon\tSplic\tNear\tFar\texpressPheAsso\n")
        for g in SNPloci.keys():
            if g == "none":
                continue
            count=SNPcount[g]
            total=count['Exonic']+count['Near']+count['Far']+count['Splic']
            of2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(g,total,count['NoSyn'],count['Exonic'],count['Splic'],count['Near'],count['Far']))

    return snp,gene
'''    

'''

def judgeCut(sign,df):
    cut = sign 
    if cut == "B":
        cut=float(1)/len(df)
    elif cut == "F":
        cut=float(0.05)/len(df)
    elif cut == "top1":
        cut = 0.5
    else:
        cut = float(cut)

    return cut


def expressPheAsso(expressDict,pheDict,gname=None,tname=None):
    x=[]
    y=[]
    for sam in pheDict.keys():
        if pheDict[sam] == "NA" or sam not in expressDict:
            continue
        else:
            x.append(float(expressDict[sam]))
            y.append(float(pheDict[sam]))
    x = sm.add_constant(x)
    model = sm.OLS(y,x).fit()

    if not gname == None:
        plt.style.use('seaborn')
        fig,ax = plt.subplots()
        ax.scatter(x, y)
        ax.set_xlabel(gname, fontsize=12)
        ax.set_ylabel(tname, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=10)
        qqline(ax, "r", x, y)
        plt.savefig('%s_vs_%s.png' %(tname,gname))

    return model


def stTest(gt,Vals):
    gts={}
    gtl=[]
    for sample in gt.keys():
        if sample not in Vals:
            continue 
        if Vals[sample] == "NA":
            continue
        gtl=accumulateList(gtl,gt[sample])
        gts=accumulateDict(gts,str(Vals[sample]),gt[sample])

    gtl=sorted(gtl)
    p1,p2,p3 = (1,1,1)
    if len(gtl) == 3:
        t1,p1 = stats.ttest_ind(np.array(listType(gts[gtl[0]],"float")),np.array(listType(gts[gtl[1]],"float")))
        t2,p2 = stats.ttest_ind(np.array(listType(gts[gtl[1]],"float")),np.array(listType(gts[gtl[2]],"float")))
        t3,p3 = stats.ttest_ind(np.array(listType(gts[gtl[0]],"float")),np.array(listType(gts[gtl[2]],"float")))
    elif len(gtl) == 2:
        t1,p1 = stats.ttest_ind(np.array(listType(gts[gtl[0]],"float")),np.array(listType(gts[gtl[1]],"float")))
    
    p=min(p1,p2,p3)
    if p1 <0.05 or  p2 <0.05 or p3 <0.05:
        return 1,p,gts
    else:
        return 0,p,gts

 
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

# Reading Function

def readSNPs(gzvcf,sitesl,out=None):
    snpsGT={}
    snpsINFO={}
    vcf_reader = VCF(gzvcf)
    samples=vcf_reader.samples
    with open("%s.sign.BasicGT.lst" %(out),"w") as of:
        of.write("Chr\tPos\tRef\tAlt\t%s\n" %("\t".join(samples)))
        for l in sitesl:
            for v in vcf_reader('%s:%s-%s' %(l[1],l[2],l[2])):
                snpsGT[l[0]] =dict(zip(samples,v.gt_types))
                snpsINFO[l[0]]=[v.REF,v.ALT[0],v.gt_phases]
                of.write("%s\t%s\t%s\t%s\t%s\n" %(l[1],l[2],v.REF,v.ALT[0],"\t".join(listType(v.gt_types.tolist(),"str"))))  #gt_types(0 ref hom,1 het, 2 unknow ,3 alt hom) 
                #of.write("%s\t%s\t%s\t%s\t%s\n" %(l[1],l[2],v.REF,v.ALT,"\t".join(v.genotypes)))
    return snpsGT,snpsINFO

def readAnnovar(vf,evf,sites):
    ginfo={}
    gcount={}
    snp2gene={}
    exon={}
    for l in open(vf,"r"):
        ls =l.strip().split()
        tags="%s#%s" %(ls[2],ls[3])
        # gene,loci
        if str(ls[2]) not in sites:
            continue
        if int(ls[3]) in sites[str(ls[2])]:
            snp2gene[tags]=[]
            pg=re.split(r'\(|,',ls[1])
            if "exonic" in ls[0]:
                gcount=accumulateDict(gcount,1,pg[0],"Exonic")
                ginfo=accumulateDict(ginfo,tags,pg[0])
                snp2gene[tags].append([pg[0],"Exonic"])
                exon["%s#%s"%(ls[2],ls[3])] = pg[0]
            elif "splicing" in ls[0]:
                gcount=accumulateDict(gcount,1,pg[0],"Splic")
                ginfo=accumulateDict(ginfo,tags,pg[0])
                snp2gene[tags].append([pg[0],"Splic"])
            elif "intergenic" in ls[0]:
                if  not pg[0] == "NONE":
                    gcount=accumulateDict(gcount,1,pg[0],"Far")
                    ginfo=accumulateDict(ginfo,tags,pg[0])
                    snp2gene[tags].append([pg[0],"Far"])
                if  not pg[1] == "NONE":
                    gcount=accumulateDict(gcount,1,pg[2],"Far")
                    ginfo=accumulateDict(ginfo,tags,pg[2])  
                    snp2gene[tags].append([pg[2],"Far"])
            else:  # intron\UTR\upstream\downstream
                if ";" in ls[1]:
                    gg=ls[1].split(";")
                    gcount=accumulateDict(gcount,1,gg[0],"Near")
                    ginfo=accumulateDict(ginfo,tags,gg[0])
                    snp2gene[tags].append([gg[0],"Near"]) 
                    gcount=accumulateDict(gcount,1,gg[1],"Near")
                    ginfo=accumulateDict(ginfo,tags,gg[1]) 
                    snp2gene[tags].append([gg[1],"Near"])
                else:   
                    gcount=accumulateDict(gcount,1,pg[0],"Near")  
                    ginfo=accumulateDict(ginfo,tags,pg[0]) 
                    snp2gene[tags].append([pg[0],"Near"])        

    for le in open(evf,"r"):
        les=le.strip().split("\t")
        tags="%s#%s"%(les[3],les[4])
        #print (tag)
        if tags in exon: 
            #print (le)
            if "nonsynonymous"  in les[1]:
                gcount=accumulateDict(gcount,1,exon[tags],"NoSyn")
                snp2gene[tags]=[[exon[tags],"NoSyn"]]

    for g in gcount.keys():
        gcount[g] = fillDictValue(gcount[g],['Near','Far','Splic','Exonic','NoSyn'],0)

    return gcount,ginfo,snp2gene

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
                    P=accumulateDict(P,v,k1,k2) 
                if order == "v":  ##  V-line  as key
                    P=accumulateDict(P,v,k2,k1)
    return P,Hlist,Vlist

def readTableDict(infile,ck=0,cv=1):
    info={}
    with open(infile,"r") as inf:
        for line in inf:
            ss = line.strip().split("\t")
            if cv != -1:
                info[ss[ck]] = ss[cv]
            else:
                info[ss[ck]] = line
    return info

# Gernal Function 
def listType(mylist,types):
    nl=[]
    for s in mylist:
        if types == "str":
            nl.append(str(s))
        elif types == "float":
            nl.append(float(s))
    return nl

def accumulateList(myList,val):
    if val not in myList:
        myList.append(val)
    return myList

def fillDictValue(myDict,keylist,val):
    for k in keylist:
        if k not in myDict:
            myDict[k] = val
    return myDict

def accumulateDict(myDict,val,keyA,keyB="no"):
    '''
        accumulate a 1~2d Dictionary
    '''
    if keyB == "no":
        if keyA in myDict:
            if val.isdigit():
                myDict[keyA] += val
            else:
                myDict[keyA].append(val)         

        else:
            if val.isdigit():
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


if __name__ == '__main__':
    logger = configure_logging()
    main()
