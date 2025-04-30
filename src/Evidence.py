#!/usr/bin/env python3 
import re
import sys 
import logging
import argparse
import numpy as np
import pandas as pd
import Common
import statsmodels.api as sm
import matplotlib.pyplot as plt
from qmplot import manhattanplot
from scipy import stats
from cyvcf2 import VCF

#logger=""

def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-i", "--input",required=True, help="Path of the Gwas or Other site-value Result of Candidata SNPs, such as *significant.tsv ",type=str)
    parser.add_argument("-p", "--phe",required=True, help="Path of the input phenotype file",type=str)
    parser.add_argument("-v", "--vcf",required=True, help='vcf file in bgzip format and has bcftools-csi index',type=str)
    parser.add_argument("-a", "--annovar",required=True,help='Prefix of snp Annovar result file, *variant_function and *exonic_variant_function',type=str)
    parser.add_argument("-f", "--function",help='Function Annotation of Genes. Support multi-column information, but the first column must be the gene name and include the table header',type=str)
    #parser.add_argument("-s", "--sign",help="Significance cut line: B = 1/n || F = 0.05/n || a self-defined float",default="B")
    parser.add_argument("-o", "--out",help="Prefix of outfiles",type=str,default="Evi")
    parser.add_argument("-w", "--weight",help='Other evidence with weight',type=str)
    parser.add_argument("--chrom",help="The column of Chromosom name in the inputfile",type=str,default="CHROM") 
    parser.add_argument("--pos",help="The column of Maker pos in the inputfile",type=str,default="POS")
    parser.add_argument("--value",help="The column of Value in the inputfile",type=str,default="P")
    parser.add_argument("--name",help="The column of SNP id in the inputfile",type=str,default="ID")
    parser.add_argument("--bed",help="Table split bed file of gene or mRNA with at least 5 column: Chrom start end strand gene_name, required if loci type file exist in --weight",type=str,default=None)
    parser.add_argument("--base_scores",help="The weights of the Basic GWAS result (range from 0~1)",type=float,default=0.8)
    parser.add_argument("--snp_scores",help="The weights of the following various type SNPs: [NoSys, Sys, Splic, Near, Far]",nargs="+",default=[10,5,3,2,1])

    parsed_args = parser.parse_args(args)
    print (parsed_args)

    Common.check_path_exists(parsed_args.input)
    Common.check_path_exists(parsed_args.phe)
    Common.check_path_exists(parsed_args.vcf)
    Common.check_path_exists(f"{parsed_args.annovar}.variant_function")
    Common.check_path_exists(f"{parsed_args.annovar}.exonic_variant_function")
    if parsed_args.function:
        Common.check_path_exists(parsed_args.function)    

    return parsed_args


def main(args=None):
    args=get_args(args)

    #logger = configure_logging()
    #logger.info('Starting main Program')
    #logger.info(args)

    out= args.out

    #logger.info(f"Reading GWAS result from {args.input}")

    df=pd.read_table(args.input, sep="\t")

    name_col= {args.chrom:"CHROM",args.pos:"POS",args.value:"P",args.name:"ID"}
    df = df.rename(columns=name_col)
    #cut=Common.judge_sign(float(1)/len(df))
    #print(df)

    #logger.info("Gwas cutline is %s"%(cut))
    (sites,sitesl,GwasP) = significanceDF(df,out)
    #print(sitesl)
    #exit(0)
    #logger.info(f"Reading Candidata SNP Genotype info from {args.vcf}")
    (cgt,cinfo)=readSNPs(args.vcf,sitesl,out)

    #logger.info(f"Reading significance SNP function from annovar result {args.annovar}.variant_function and {args.annovar}.exonic_variant_function") 
    (SNPcount,SNPloci,snp2gene) = readAnnovar(f"{args.annovar}.variant_function",f"{args.annovar}.exonic_variant_function",sites) 

    #logger.info("Reading Phenotype data from %s" %(args.phe)) 
    Phe = readTableDict(args.phe,0,2)

    (snp,candidate_genes)=outBasicModel(out,sites,GwasP,cgt,cinfo,SNPcount,SNPloci,snp2gene,Phe,args.function) 

    
    if args.weight:  # with eight model
        main_df = pd.read_table(f"{out}.sign.BasicGene.xls",index_col="ID")
        df_handle = main_df.iloc[:, :8].copy()
        df_function = main_df.iloc[:, 8:]
        #print(df_handle)
        #print(df_function)

        evidens=Common.read_file(args.weight,mode="list",vals=[0,1,2,3])
        evidens = sorted(evidens, key=lambda x: 0 if x[1] == 'eqtl' else 1)

        second_score_dict={}
        print(evidens)
        for name,model,value,file in evidens:
            if model == "express":
                (df_handlese)=evidence_express(file,df_handle,name,Phe)     # Phe vs FPKM
                #second_score_dict.append([name,value])    
            elif model == "gene":   # Known Candidate Gene
                df_handle=evidence_genes(file,df_handle,name)
                #second_score_dict.append([name,value])          
            elif model == "eqtl":   # Known SNP-Gene associate
                df_handle=evidence_eqtl(file,df_handle,name,basic_site_file=f"{out}.sign.BasicSite.xls")
                #second_score_list.append([f"{name}_global",value])               
                #second_score_dict[name]=name,value
            elif model == "loci":   # Known Candidate Gene
                df_handle=evidence_loci(file,df_handle,name,args.bed)
                #second_score_dict[name]=value
            else:
                print(f"# Unknow evidence data type {model},Ignore and skip file {file}")
            second_score_dict[name]=value

        #df_handle.to_csv(f"{out}.sign.Weight.genes.txt",sep="\t")
        print(df_handle)
        print(second_score_dict)


        base_score_dict = dict(zip(["NoSys","Sys","Splic","Near","Far"],args.snp_scores))
        df_handle = cal_df_score(df_handle,base_score_dict,"Base_score")

        df_handle = cal_df_score(df_handle,second_score_dict,"Evi_score")

        add_dict=dict(zip(["Base_score_normalized","Evi_score_normalized"],[args.base_scores,1-args.base_scores]))
        print(add_dict)
        df_handle = cal_df_score(df_handle,add_dict,"Total_score")

        df_handle.sort_values(by='Total_score_normalized', ascending=False, inplace=True)
        print(df_handle)

        if args.function:
            df_handle=add_function(df_handle,args.function)
        print(df_handle)
        df_handle.to_csv(f"{out}.sign.Weight.genes.xls",sep="\t")


def add_function(df,function_file):
    functions,header=Common.read_file(function_file,mode="dict",vals=[],keys=[0],header=True,sep="\t")
    for i in range(1,len(header)):
        info_dict={}
        for g in functions.keys():
            info_dict[g] = functions[g][i]
        df=updata_df(df,header[i],info_dict,missing="---")

    return df


def cal_df_score(df,weights,name):
    weight_dict = {col: float(val) for col, val in weights.items()}
    df[f'{name}'] = sum(df[col] * weight for col, weight in weight_dict.items())
    df[f'{name}_normalized'] = df[name] / df[name].sum() *100
    df[f'{name}_normalized'].round(4)
    return df


def evidence_loci(file,df_handle,name,bedfile):
    if not bedfile:
        print("# Error: Please specify the bed file containing mRNA or genes using the --bed !")
        exit(1)

    gene_location = Common.read_file(bedfile,mode="dict",keys=[4],vals=[0,1,2,3])
    #print(gene_loci)
    #locis = Common.read_file(file,mode="dict",vals=[0,1,2,3],keys=[0])   # chrom start end id
    locis = Common.read_file_accumulateDict(file,vals=[0,1,2,3],key1=0)

    gene_loci={}
    candidate_genes=df_handle.index.tolist()
    for gene in candidate_genes:
        #print(gene)
        chrid=gene_location[gene][0]
        overlapping_regions=Common.search_cloesd_region(gene_location[gene],locis[chrid])
        if len(overlapping_regions) >= 1:
            gene_loci[gene] = 1
        else:
            gene_loci[gene] = 0

        #print(gene,overlapping_regions)
    df_handle=updata_df(df_handle,name,gene_loci)
    return df_handle


def evidence_eqtl(file,df_handle,name,basic_site_file="out.sign.BasicSite.xls"):
    basic_snps = Common.read_file_accumulateDict(basic_site_file,vals=[0],key1=3)
    eqlt_snps  = Common.read_file_accumulateDict(file,vals=[0],key1=1)

    total_snps = [x[0] for x in Common.read_file(basic_site_file,mode="list",vals=[0])]
    snps_eqlt  = Common.read_file_accumulateDict(file,vals=[1],key1=0)

    for snp in total_snps:
        if not snp in snps_eqlt:
            continue
        for gene in snps_eqlt[snp]:
            if gene not in df_handle.index:
                basic_snps[gene]=[]
                print(f"# Add new candidate_genes {gene} according eqtl file {file} ")
                df_handle.loc[gene] = {
                    'TotalSNP': 0,
                    'NoSys': 0,
                    'Sys': 0,
                    'Splic': 0,
                    'Near': 0,
                    'Far': 0,
                    'PheGT-test-number': 0,
                    'PheGT-test-ratio': 0,
    }

    gene_snps_local={}
    gene_snps_global={}
    candidate_genes=df_handle.index.tolist()
    #print(candidate_genes)
    for gene in candidate_genes:
        common_snps = []
        comon_all_snps = []
        if gene in eqlt_snps:
            common_snps = set(basic_snps[gene]) & set(eqlt_snps[gene]) #  local-eqtl
            comon_all_snps = set(total_snps) & set(eqlt_snps[gene])
        gene_snps_local[gene] = len(common_snps)
        gene_snps_global[gene] = len(comon_all_snps)

    df_handle=updata_df(df_handle,name,gene_snps_global)
    #df_handle=updata_df(df_handle,f"{name}_local",gene_snps_local)
    #df_handle=updata_df(df_handle,f"{name}_global",gene_snps_global)

    return df_handle


def evidence_genes(file,df_handle,name):
    gene_exist={}
    genes_tmp = Common.read_file(file,mode="list",sep="\t")
    genes = [x[0] for x in genes_tmp]
    
    candidate_genes=df_handle.index.tolist()
    for gene in candidate_genes:
        if gene in genes:
            gene_exist[gene]=1
        else:
            gene_exist[gene]=0

    df_handle=updata_df(df_handle,name,gene_exist)
    return df_handle


def evidence_express(file,df_handle,name,Phe):
    #logger.info("Reading the gene expression matrix from %s" %(args.express))
    gene_asso={}
    G,Gsample,Genes  = readData(file,"v")
    candidate_genes = df_handle.index.tolist()
    for gene in candidate_genes:
        if gene not in G:
            print("%s has no Transcript Data, skipping" % (gene))
            gene_asso[gene] = "NA"
        else:
            lineregress=expressPheAsso(G[gene],Phe)
            #print(gene,lineregress.pvalues,lineregress.pvalues[1])
            p = 0
            if not np.isnan(lineregress.pvalues[1]) and lineregress.pvalues[1] <0.01:
                p = 1
            #print(gene,lineregress.pvalues,lineregress.pvalues[1],p)
            gene_asso[gene] = p

    df_handle=updata_df(df_handle,name,gene_asso)

    return df_handle


def updata_df(df,index_name,data_dict,missing="None"):
    df[index_name] = df.index.map(data_dict)
    df[index_name].fillna(missing, inplace=True)
    #new_order = [c for c in df.columns if c != 'B'] + ['B']
    return df



def outBasicModel(out,sites,GwasP,cgt,cinfo,SNPcount,SNPloci,snp2gene,Phe,function_file):
    snp=[]
    gene=[]
    sign_snp_number={}
    #logger.info("Outputing significance SNP function, GWAS-Pvalue")

    # needed snp2gene sites cinfo GwasP cgt
    with open("%s.sign.BasicSite.xls" % (out),"w") as of1:
        of1.write("#ID\tRef\tAlt\tType\tGene\tP-GWAS\tP-CorPhe\n")
        for snp in snp2gene.keys():
            rs = snp.split("#")
            if rs[0] not in sites or int(rs[1]) not in sites[rs[0]]:
                #logger.warning("Check SNP %s %s" %(rs[0],rs[1]))
                continue
            sID = sites[rs[0]][int(rs[1])]
            (pp1,p1,r)=stTest(cgt[sID],Phe)     # student t-test SNP vs Phe
            for recode in snp2gene[snp]:
                if p1 < 0.01 :
                    sign_snp_number=accumulateDict(sign_snp_number,1,recode[0])
                of1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(sID,cinfo[sID][0],cinfo[sID][1][0],recode[0],recode[1],GwasP[sID],p1))

    #logger.info("Outputing Genes-SNP statistical table")

    if function_file:
        functions,function_header = Common.read_file(function_file,mode="dict",vals=[],keys=[0],sep="\t",header=True)
        print(function_header) 
    
    out_list=[]

    out_header=["ID","TotalSNP","NoSys","Sys","Splic","Near","Far","PheGT-test-number","PheGT-test-ratio"]
    if function_file:
        out_header += function_header[1:]

    for g in SNPloci.keys():
        if g == "NONE":
            continue

        count=SNPcount[g]
        total=count['Exonic']+count['Near']+count['Far']+count['Splic']

        if g not in sign_snp_number:
            sign_snp_number[g]=0
        sign_snp_rate = round(sign_snp_number[g]/total,2)

        line = [g,total,count['NoSyn'],count['Exonic']-count['NoSyn'],count['Splic'],count['Near'],count['Far'],sign_snp_number[g],sign_snp_rate]

        if function_file:
            if g not in functions:
                functions[g]=["---" for x in function_header]
            line += functions[g][1:]

        out_list.append(line)

    out_list=sorted(out_list,key=lambda x: x[1],reverse=True)
    Common.write_list_to_file(out_list,f"{out}.sign.BasicGene.xls",header=out_header,vals=[],sep="\t")
    #with open("%s.sign.BasicGene.xls" % (out),"w") as of2:

    return snp,[x[0] for x in out_list]
'''    

'''

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

 
def significanceDF(df,out):
    sites={}
    GwasP={}
    sitesl=[]

    #signDf=df[df['P'] < cut ]
    signDf=df
    # out signDf  
    #logger.info(f"{len(signDf)} significance SNP were founded !")
    #signDf.to_csv("%s.sign.sites" %(out),index=False,sep="\t")
 
    for index,row in signDf.iterrows():
        c=str(row['CHROM'])
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
                if  not pg[2] == "NONE":
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



def configure_logging():
    #logging.basicConfig(filename="EGCV.log",level=logging.INFO,format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",)
    return logging.getLogger("EVCG")


if __name__ == '__main__':
    #logger = configure_logging()
    main()
