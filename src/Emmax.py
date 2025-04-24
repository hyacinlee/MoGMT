#!/usr/bin/env python3
import Common
import Visualize
import ParaRun
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse


def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-v", "--vcf",required=True, help="Path of the input vcf file")
    parser.add_argument("-p", "--phe",required=True, help="Path of the input phenotype file")
    parser.add_argument("-c","--cov",help="Covariance, could be PCA, structure or other files",type=str) 
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out") 
    parser.add_argument("-m","--maf",help="Exclude variants with minor allele frequency lower than this",type=float,default="0.05")
    parser.add_argument("-s","--sign",help="Significance cut line: B = 1/n || F = 0.05/n || a self-defined float",default="B")
    parser.add_argument("-e", "--emmax",help="Path of the emmax dir,which include emmax-kin-intel64 and emmax-intel64 ",default="/home/minghui/software/reseq/EMMAX/")
    parser.add_argument("-l", "--plink",help="Path of the plink v1.9",default="/home/minghui/software/reseq/plink1.9/plink")
     
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args
  
def main(args=None):
    args=get_args(args)

    Common.check_path_exists(f"{args.emmax}/emmax-intel64")
    Common.check_path_exists(f"{args.emmax}/emmax-kin-intel64")
    Common.check_path_exists(f"{args.plink}")

    samples=process_vcf_phe(args)

    if args.cov:
        process_covariance(args.cov,samples)

    traits=process_gwas(args)

    #snps = pd.read_csv("GWAS.map",header=None,sep="\t",names=["chrom","snp_id","unknow","pos"])
    #snps.set_index("snp_id", inplace=True)                 # panda df速度更慢
    #snps=Common.read_file("GWAS.map","dict",vals=[0,3],keys=[1])    
    #significant=Common.judge_significant(sign,snps)
    snps=Common.read_file("GWAS.snpID","dict",vals=[0,1],keys=[2])    
    significant=Common.judge_significant(args.sign,snps)
    with open("GWAS.threshold","w") as f:
        f.write(f"{significant}\n")
 
    if not Path("GWAS.outputs.done").exists():
        process_outputs(traits,snps,significant)

    if not Path("GWAS.mahantun.done").exists():
        process_manhantun(traits,significant)


def process_manhantun(traits,significant):
    for trait in traits:
        Visualize.ManhantanMain(f"Traits.{trait}.emmax.ps.tsv",out=f"Traits.{trait}.emmax.ps.manhantun",sign=significant)
    Common.run_command(f"touch GWAS.manhantun.done")

def process_outputs(traits,snps,significant):
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(process_trait_output, trait, snps, significant) for trait in traits]

    for future in as_completed(futures):
        future.result() 

    Common.run_command(f"touch GWAS.outputs.done")

    #return significant


def process_trait_output(trait,snps,significant):
    with open(f"Traits.{trait}.emmax.ps", "r") as f,open(f"Traits.{trait}.emmax.ps.tsv", "w") as f1,open(f"Traits.{trait}.emmax.ps.significant.tsv", "w") as f2:
        f1.write("CHROM\tPOS\tID\tP\n")
        f2.write("CHROM\tPOS\tID\tP\n")
        #data=[]
        for l in f:
            ls = l.strip().split("\t")
            #chrom, pos = snps.loc[ls[0], ["chrom", "pos"]]
            chrom, pos = snps[ls[0]]
            f1.write(f"{chrom}\t{pos}\t{ls[0]}\t{ls[3]}\n")
            if float(ls[3]) <= significant:
                f2.write(f"{chrom}\t{pos}\t{ls[0]}\t{ls[3]}\n")

    #Manhantun.ManhantanMain(f"Traits.{trait}.emmax.ps.tsv",out=f"Traits.{trait}.emmax.ps.manhantun",sign=significant)

    #Manhantun.manhatan_fig(df,significant,ylab="-log10(P)",out="Traits.{trait}.manhatun",sub=None,hightlight=None)
    #
        # use panda df 
            #row={'#CHROM': chrom,'POS': int(pos),'ID': str(ls[0]),'P':float(ls[3])}
            #data.append(row)
        #df = pd.DataFrame(data)            
        #signDf=df[df['P'] < significant ]
        #df.to_csv(f"Traits.{trait}.emmax.ps.tsv",index=False,sep="\t")
        #signDf.to_csv(f"Traits.{trait}.emmax.ps.significant.tsv",index=False,sep="\t")



def process_gwas(args):
    data,traits,samples = Common.read_matrix_data(f"GWAS.phe","h")
    print (f"Run GWAS with population size: {len(samples)} with {len(traits)} Traits")
    print (f"Traits for run: {traits}")

    cmds=[]
    for trait in traits:
        trait_phes,samples = Common.sort_dict_by_list(data[trait],samples)
        Common.write_list_to_file(zip(samples,samples,trait_phes),f"Traits.{trait}.phe",vals=[])
        cmds.append(f"{args.emmax}/emmax-intel64 -v -d 10 -t GWAS -k GWAS.aBN.kinf -p Traits.{trait}.phe -o Traits.{trait}.emmax -c GWAS.covariance")
    
    Common.write_list_to_file(cmds,"./cmd.GWAS.sh")
    ParaRun.runlocal("./cmd.GWAS.sh",lines=1,threads=8)

    return traits


def process_vcf_phe(args):
    if not Path("GWAS.prepared.done").exists():
        with open(args.vcf,"r") as f1, open("GWAS.snpID","w") as fo:
            for line in f1:
                ls = line.strip().split("\t")
                if line.startswith("#CHR"):
                    vcf_samples = line.strip().split("\t")[9:]
                elif not line.startswith("#")
                    if not ls[2] == "."
                        fo.write(f"{ls[0]}\t{ls[1]}\t{ls[2]}\n")
                    else:
                        print(f"# There is no SNP-ID infomation cols in {args.vcf}, you can use the Vcftools module to add it")
                        exit(1)

        (phe_data,headinfo)=Common.read_file(args.phe,mode="dict",keys=[0],header=True)  
        (phe_data,samples)=Common.sort_dict_by_list(phe_data,vcf_samples)

        Common.write_list_to_file(phe_data,"GWAS.sample",vals=[0,0])
        Common.write_list_to_file(phe_data,"GWAS.phe",vals=[],header=headinfo)

        Common.run_command(f"{args.plink} --vcf {args.vcf} --allow-extra-chr --make-bed --out GWAS --keep GWAS.sample --maf {args.maf}")
        Common.run_command(f"{args.plink} --bfile GWAS --out GWAS --output-missing-genotype 0 --recode 12 transpose")
        Common.run_command(f"{args.emmax}/emmax-kin-intel64 -v -d 10 GWAS")   # Generate kinship matrix
        Common.run_command(f"touch GWAS.prepared.done")

    samples=Common.read_file("GWAS.sample","list",vals=[0])
    return samples


def process_covariance(cov_file,samples):
    """Process structure/PCA file"""
    cov_data=Common.read_file(cov_file,mode="list",vals=[]) # read cov without header
    cov_data_new=[] 
    for ls in cov_data:
        ls.insert(1,1)
        ls.insert(1,ls[0])
        cov_data_new.append(ls)
    Common.write_list_to_file(cov_data_new,"GWAS.covariance",vals=[])



if __name__ == "__main__":
    main()