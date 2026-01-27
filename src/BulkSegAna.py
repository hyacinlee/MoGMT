#!/usr/bin/env python3
import os
import sys
import Common
import Visualize
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse

def get_args(args):
    parser = argparse.ArgumentParser(
    formatter_class = Common.CustomFormatter,
    description = '''
RunCmd :
    Funciton:  Performing Bulk-sequencing Analysis using Euclidean-distance or SNP-index.
    Writer:    Meng Minghui 
'''
    )
    parser.add_argument('-v','--vcf',help='vcf file',required=True,type=str)
    parser.add_argument('-t','--type',help='Use Euclidean-distance or SNP-index',type=str,choices=["ed","index"],required=True)
    parser.add_argument('-o','--out',help='pre-name of out files',type=str,default="BSA")
    parser.add_argument('-s','--sign',help="Significance cut line: B = 1/n || F = 0.05/n || top1 = top 0.01 || top5 = top 0.05 || a self-defined float",nargs="+",default=["top1"])
    parser.add_argument('-step',help='Step-size of sliding window, kb',type=int,default=None)
    parser.add_argument('-wind',help='Wind-size of sliding window , kb',type=int,default=None)
    parser.add_argument('-fit',help='Fit method of sliding window',choices=["mean","linear"],type=str,default="mean")
    parser.add_argument('-minDepth',help='min Depht of each bulks',type=int,default=10)
    parser.add_argument('-minSite',help='min site of a window',type=int,default=10)
    parser.add_argument('-s1',help='Name of Segregation Bulk1',type=str,default=None,required=True)
    parser.add_argument('-s2',help='Name of Segregation Bulk2',type=str,default=None,required=True)
    parser.add_argument('-p1',help='Name of Parent Bulk1 ',type=str,default=None)
    parser.add_argument('-p2',help='Name of Parent Bulk2',type=str,default=None)
    
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args 

def main(args=None):
    args=get_args(args)
    main_BSA(args)


def main_BSA(args):

    site_out= f"{args.out}.site.txt"

    if not os.path.exists(site_out):
        out1=open(site_out,"w")
        if args.type == "ed":
            out1.write("Chrom\tPos\tID\tRef\tAlt\tValue\ted4\tBulk1\tBulk2\n")
        elif args.type == "index":
            out1.write("Chrom\tPos\tID\tRef\tAlt\tValue\tIndex1\tIndex2\tBulk1\tBulk2\n")
        name_index = {}

        print("# Start cal site BSA.... ")
        for line in open(args.vcf):
            if line.startswith("##"):
                continue
            elif line.startswith("#CHR"):
                for i, name in enumerate(line.strip().split()):
                    name_index[name] = i  
            else:
                info = line.strip().split()
                if args.type == "ed":
                    result = calED(info,name_index[args.s1],name_index[args.s2],args.minDepth)
                elif  args.type == "index":
                    result = calIndex(info,name_index[args.s1],name_index[args.s2],args.minDepth)
                    
                if len(result) == 0 :
                    continue
                out1.write("\t".join(map(str,result)) + "\n") 

    else:
        print("# result file {site_out} exist, remove it if you want re-run site-BSA")

    Site = pd.read_table(f"{args.out}.site.txt", sep="\t")
    #df_plot = df.rename(columns={"ed2":"Value"}) 
    cutoff=Common.judge_significant(args.sign,Site,"Value")
    Visualize.manhatan_fig_v1(Site,cutoff,"ED2",f"{args.out}.site")
    output_candidate(Site,"Value",f"{args.out}.site",args.sign,cutoff)

   
    if args.wind:
        Wind=sliding_window_fit(Site,f"{args.out}.wind{args.wind}k.txt",args.step*1000,args.wind*1000,args.minSite,method=args.fit,chrom_col='Chrom', pos_col='Pos', value_col='Value')   
        cutoff1=Common.judge_significant(args.sign,Wind,"Value")
        draw_combined(Site,Wind.fillna(0),cutoff1,f"{args.out}.wind{args.wind}k")
        candidate=output_candidate(Wind,"Value",f"{args.out}.wind{args.wind}k",args.sign,cutoff1,mergeOut=True)


def output_candidate(df,colum,out,sign,cuts,mergeOut=False):
    for i,c in enumerate(cuts):
        df_c = df[df[colum] > cuts[i]]
        df_c.to_csv(f"{out}.candidate.{sign[i]}.txt",index=False,sep="\t")
        print(f"# Output candidate sites/regions to {out}.candidate.{sign[i]}.txt ")
        if mergeOut:
            merged_candidate=Common.merge_overlapping_intervals(df_c)
            merged_candidate.to_csv(f"{out}.candidate.{sign[i]}.merged.txt",index=False,sep="\t")


def draw_combined(Site,Wind,cuths,outfile,posC="Value",lineC="Value"):
    xmax=Site['Pos'].max()
    ymax=Site[posC].max()

    bins=list(Wind['Chrom'].unique())
    fig,axes = plt.subplots(len(bins),1,sharex=True,figsize=(2*len(bins),len(bins)))
    for i in range(len(bins)):
        name = bins[i]
        axes[i].text(xmax*1.05,ymax/2,name,fontsize=15,verticalalignment="top",horizontalalignment="left")
        axes[i].scatter(Site[Site['Chrom']==name]['Pos'],Site[Site['Chrom']==name][posC],s=0.5,c='lightgray')
        axes[i].plot(Wind[Wind['Chrom']==name]['Start'],Wind[Wind['Chrom']==name][lineC],c="firebrick")
        axes[i].spines['top'].set_color("white")
        axes[i].spines['right'].set_color("white")
        axes[i].xaxis.set_major_formatter(plt.FuncFormatter(millions))

        for c in cuths:
            axes[i].axhline(y=c,ls="--",c="navy",lw=1)

    plt.savefig(fname="%s.png" % (outfile),dpi=300)
    plt.show()



def sliding_window_fit(df,out,step_size,wind_size,minSite,method="mean",chrom_col='Chrom', pos_col='Pos', value_col='Value'):

    results = []
    out=open(out,"w")
    out.write(f"Chrom\tStart\tEnd\tValue\tNumber\n")
    total_winds=0
    vaild_winds=0

    for chrom in df[chrom_col].unique():
        chrom_data = df[df[chrom_col] == chrom].copy()
        chrom_data = chrom_data.sort_values(pos_col).reset_index(drop=True)
        
        if len(chrom_data) == 0:
            continue
            
        positions = chrom_data[pos_col].values
        values = chrom_data[value_col].values
        
        min_pos = positions.min()
        max_pos = positions.max()
        
        left_ptr = 0  
        window_id = 0
        window_start = 0
        
        while window_start <= max_pos:
            total_winds += 1
            window_end = window_start + wind_size
            
            while left_ptr < len(positions) and positions[left_ptr] < window_start:
                left_ptr += 1
            
            right_ptr = left_ptr
            while right_ptr < len(positions) and positions[right_ptr] < window_end:
                right_ptr += 1
            
            snp_indices = slice(left_ptr, right_ptr)
            snp_number = right_ptr - left_ptr
            
            fitted_value = None
            if snp_number >= minSite:
                vaild_winds += 1 
                fitted_value = fit_values(positions[snp_indices],values[snp_indices],method)
                #window_positions = positions[snp_indices]
                #window_values = values[snp_indices]
                    
            out.write(f"{chrom}\t{window_start}\t{window_end}\t{fitted_value}\t{snp_number}\n")
            #print(chrom,window_start,window_end,fitted_value,snp_number)
            results.append({
                'Chrom': chrom,
                'Start': window_start,
                'End': window_end,
                'Value': fitted_value,
                'Number': snp_number
            })

            
            window_start += step_size
            window_id += 1
    print(f"# Get {vaild_winds} vailded windows out of {total_winds}")

    return pd.DataFrame(results)
    

def fit_values(window_positions,window_values,method):
    if method == "mean":
        fitted_value = np.mean(window_values)
    
    elif method == "linear":
        if snp_number > 1:
            slope, intercept, _, _, _ = stats.linregress(window_positions, window_values)
            mid_pos = (window_start + window_end) / 2
            fitted_value = intercept + slope * mid_pos
        else:
            fitted_value = window_values[0]
    
    else:
        raise ValueError("method must be 'mean' or 'linear'")

    return fitted_value



def calED(info,s1,s2,minDepth):
    ds1 = info[s1].split(":")[1]
    ds2 = info[s2].split(":")[1]

    a1,a2 = map(int,ds1.split(","))
    b1,b2 = map(int,ds2.split(","))
    total =  a1 + a2 + b1 + b2
    if a1+a2 < minDepth or b1+b2 < minDepth:
        return []
    else:
        ed2 = (float(a1)/(a1+a2)-float(b1)/(b1+b2))**2 + (float(a2)/(a1+a2)-float(b2)/(b1+b2))**2
        ed4 = ed2**2
        res= [info[0],info[1],info[2],info[3],info[4],ed2,ed4,ds1,ds2]
        return res


def calIndex(info,s1,s2,minDepth):
    ds1 = info[s1].split(":")[1]
    ds2 = info[s2].split(":")[1]

    a1,a2 = map(int,ds1.split(","))
    b1,b2 = map(int,ds2.split(","))
    total =  a1 + a2 + b1 + b2
    if a1+a2 < minDepth or b1+b2 < minDepth:
        return []
    else:
        ed2 = (float(a1)/(a1+a2)-float(b1)/(b1+b2))**2 + (float(a2)/(a1+a2)-float(b2)/(b1+b2))**2
        ed4 = ed2**2
        res= [info[0],info[1],info[2],info[3],info[4],ed2,ed4,ds1,ds2]
        return res


def millions(x,pos):
    return '{:1.1f}M'.format(x*1e-6)


if __name__ == '__main__':
    main()
