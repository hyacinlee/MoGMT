#!/usr/bin/env python3 
import sys
import Common
import argparse
import numpy as np
import pandas as pd



def get_args(args):
    parser = argparse.ArgumentParser(description="Co localization with mutiple traits.",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-i","--input", required=True, help="Traits list file to proferm Co localization, one trait per line")
    parser.add_argument("-r","--rename", help="Re name the traits, -i must have two columns",action='store_true') 
    parser.add_argument("-d","--dir", help="Path contain *BasicGene.xls",default="./")
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="Co_localization") 
    parsed_args = parser.parse_args(args)
    
    return parsed_args
    
#corr_threshold=0.8, pc_threshold=0.8

def main(args=None):

    args=get_args(args)

    traits={}

    if args.rename:
        traits=Common.read_file(args.input,mode="dict",keys=[0],vals=[1]) 
    else:
        traits=Common.read_file(args.input,mode="dict",keys=[0],vals=[0])

    #print(traits)
    df_list = []
    for t in traits.keys():
        tf = f"{args.dir}/{t}.sign.BasicGene.xls"
        
        if not Common.check_path_exists(tf,"w"):
            continue

        dt=pd.read_csv(tf,sep="\t")
        dt.insert(0, 'Triat', t)
        if args.rename:
            dt.insert(1, 'Triat_annot', traits[t])
        df_list.append(dt)

    df = pd.concat(df_list, ignore_index=True)
    #print(df)
    df.to_csv(f"{args.out}.co_localization.data.tsv",sep="\t",index=False)  

    col_start = df.columns.get_loc("Nearest_SNP_P") + 1
    func_cols = df.columns[col_start:]

    gene_info = df.groupby("Gene")[list(func_cols)].first()

    trait_summary = df.groupby("Gene").apply(summarize_traits)

    result = (
        trait_summary
        .join(gene_info)  # 同为 Gene 索引
        .reset_index()
        .sort_values("Triat_Count", ascending=False)
    )
    #print(result)
    result.to_csv(f"{args.out}.co_localization.merge.tsv",sep="\t",index=False)  



def summarize_traits(group):
    triats = sorted(set(group["Triat"]))
    result = {
        "Triat_Count": len(triats),
        "Triats": ",".join(triats),
    }

    if "Triat_annot" in group.columns:
        # 构建 Triat -> Triat_annot 映射
        triat_map = dict(zip(group["Triat"], group["Triat_annot"]))
        annot_list = [triat_map[t] for t in triats]
        result["Triat_annots"] = ",".join(annot_list)

    return pd.Series(result)



if __name__ == "__main__":
    main()
