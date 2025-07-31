#!/usr/bin/env python3 
import sys
import Common
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer


def get_args(args):
    parser = argparse.ArgumentParser(description="Processing phenotypic data, including clustering, dimensionality reduction, etc.",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-p","--phe", required=True, help="Path of the input phenotype file")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-g","--globals",help="Global PCA and automatically cluster dimensionality reduction",action='store_true')
    group.add_argument("-l","--locals",help="Local PCA and automatically cluster dimensionality reduction with specified traits")
    parser.add_argument("-i","--impute",help="Impute the missing value use KNNImputer method",action='store_true')
    parser.add_argument("-n","--n_neighbors",help="The number of neighbors to run KNNImputer,default=5",default=5,type=int)
    parser.add_argument("-c","--corr_threshold",help="Threshold of the cor bewteen two traits to divide into a same group,default=0.6",default=0.6,type=float)
    parser.add_argument("-z","--pc_threshold",help="Threshold of the PC1 bewteen two traits to divide into a same group,default=0.8",default=0.8,type=float)
    parser.add_argument("-m","--miss_threshold",help="The maximum allowable proportion of missing values. Those exceeding this proportion will be discarded,default=0.2",type=float,default=0.2)
    parser.add_argument("-f","--fig",help="Also output the groups corr-heatmap",action='store_true')
    parser.add_argument("--subTraits",help="Extract the traits in the file.",default=None,type=str)
    parser.add_argument("--subSamples",help="Extract the Samples in the file.",default=None,type=str)
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out")   
    parsed_args = parser.parse_args(args)
    

    return parsed_args
    


#corr_threshold=0.8, pc_threshold=0.8

def main(args=None):

    args=get_args(args)

    phe_data=load_data(args.phe)
    #print(phe_data)

    if args.subTraits or args.subSamples:
        print("Get Sub df with samples or traits")
        df_sub=sub_sample_traits(phe_data,args.subTraits,args.subSamples)
        df_sub.to_csv(f"{args.out}.sub_set.tsv",sep="\t")
        correlation_matrix = df_sub.corr()
        correlation_matrix.to_csv(f"{args.out}.sub_set.cor_matrix.tsv",sep="\t")

        exit(0)


    df_mis = missing_values(phe_data, args.miss_threshold,args.out)

    df_mis_inpute = df_mis
    if args.impute:
        df_mis_inpute=impute_values(df_mis,args.n_neighbors,args.out)

    if args.globals:
        groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold)
        output(df_mis_inpute,groups,pc_values,group_dfs,group_pc_dfs,args.out,args.fig)
        cor_fig(df_mis_inpute,f"{args.out}.sample")
        cor_fig(df_mis_inpute.T,f"{args.out}.trait")
    else: #locals
        start_list=pd.read_csv(args.locals, sep='\t').iloc[:, 0].tolist()
        groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold,start_list)
        output(df_mis_inpute,groups,pc_values,group_dfs,group_pc_dfs,args.out,args.fig)




def sub_sample_traits(df,subTraits,subSamples):
    #print(df.columns)
    if subTraits:
        trait_list=Common.read_file(subTraits,mode="list",vals=[0])
        print(trait_list)
        missing_traits = [t for t in trait_list if t not in df.columns]
        sub_traits = [t for t in trait_list if t in df.columns]
        if missing_traits:
            print(f"The following traits were not found in the DataFrame and will be ignored: {missing_traits}")
        df = df.loc[:,sub_traits]
    if subSamples:
        samples_list=Common.read_file(subSamples,mode="list",sep=None)
        #print(samples_list)
        df = df.loc[samples_list]
    return df


# 进行性状分组
def group_traits(df, corr_threshold=0.8, pc_threshold=0.8, start_list=None):
    corr_matrix = df.corr()  # 保留正负号
    all_traits = set(df.columns)
    remaining_traits = all_traits.copy()

    grouped_traits = []
    group_pc_values = []
    group_dfs = []
    group_pc_dfs = []

    start_queue = list(start_list) if start_list else []

    while remaining_traits:
        if start_list:   # local
            if start_queue:
                base_trait = start_queue.pop()
                if base_trait not in remaining_traits:
                    continue
                print (f"#Cluster starts with {base_trait}")
            else:
                break
        else: # Global
            base_trait = remaining_traits.pop()

        group = [base_trait]

        if base_trait in remaining_traits:
            remaining_traits.remove(base_trait)

        #只考虑正相关性
        candidates = {
            trait for trait in remaining_traits
            if corr_matrix.loc[base_trait, trait] > corr_threshold
        }

        # 按 group 中已有性状的平均正相关性降序排序
        candidates = sorted(
            candidates,
            key=lambda trait: np.mean([corr_matrix.loc[trait, t] for t in group]),
            reverse=True
        )

        for trait in candidates:
            temp_group = group + [trait]
            temp_df = df[temp_group].dropna()

            if temp_df.shape[0] > 1:
                pc1_ratio, _ = compute_pca(temp_df)
                if pc1_ratio >= pc_threshold:
                    group.append(trait)
                    remaining_traits.remove(trait)

        group_df = df[group].dropna()
        if group_df.shape[1] == 0:
            continue

        pc1_ratio, pc1_df = compute_pca(group_df)
        grouped_traits.append(group)
        group_pc_values.append(pc1_ratio)
        group_dfs.append(group_df)
        group_pc_dfs.append(pc1_df)

    sorted_groups = sorted(
        zip(grouped_traits, group_pc_values, group_dfs, group_pc_dfs),
        key=lambda x: len(x[0]),
        reverse=True
    )

    return zip(*sorted_groups)

def output(df,groups,pc_values,group_dfs,group_pc_dfs,outprefix,fig):
    out_path  = f"{outprefix}.Combine.Triat.txt"
    out=open(out_path,"w")

    df1 = pd.DataFrame(index=df.index)
    for i, pc1_scores in enumerate(group_pc_dfs, start=1):
        df1[f"Group{i}"] = pc1_scores["PC1"]    

    # 输出分组结果
    for i, (group, pc1,dfg,dfpc) in enumerate(zip(groups,pc_values,group_dfs,group_pc_dfs), start=1):
        if len(group) > 2 and fig:
            cor_fig(dfg,f"{outprefix}.Combine.Group%s.png"%(i))     
        #print(f"Group{i}: {group}, PC1 贡献率: {pc1:.2%}")       
        gg=",".join(group) 
        out.write(f"Group{i}\t{len(group)}\t{pc1:.2%}\t{group[0]}\t{gg}\n")     # mian out 


    selected_traits = [trait for group in groups if len(group) >= 2 for trait in group]
    df3 = df[selected_traits]

    selected_traits_df4 = [group[0] for group in groups if len(group) == 1]
    df4 = df[selected_traits_df4]

    df1.to_csv(f"{outprefix}.PC1_scores.txt",sep="\t", na_rep="NA")  # 每个样本对应每个group的PC1值
    df3.to_csv(f"{outprefix}.Trait.inGroup.txt",sep="\t", na_rep="NA")  # 仅包含 group 长度 >= 2 的性状
    df4.to_csv(f"{outprefix}.Trait.unGroup.txt",sep="\t", na_rep="NA")


def load_data(file_path):   
    df = pd.read_csv(file_path,sep="\t", index_col=0)  # 以第一列（样本名）作为索引
    df.replace([np.inf, -np.inf], np.nan, inplace=True)  # 处理 inf
    return df

def missing_values(df,miss_threshold,out):
    print(f"The shape of input data: {df.shape}")
    df = df.select_dtypes(include=[np.number]) 
    missing_rates = df.isna().mean()
    columns_to_drop = missing_rates[missing_rates > miss_threshold].index
    df.drop(columns=columns_to_drop, inplace=True) 
    print(f"Deletes {len(list(columns_to_drop))} traits column with a deletion rate greater than {miss_threshold*100:.0f}% : {list(columns_to_drop)}")
    df.to_csv(f"{out}.dropmis.txt",sep="\t", na_rep="NA")
    return df

def impute_values(df,n_neighbors,out):
    print(f"The shape of input data: {df.shape}")
    imputer = KNNImputer(n_neighbors=n_neighbors)
    df_imputed = pd.DataFrame(imputer.fit_transform(df), columns=df.columns, index=df.index)
    df_imputed.to_csv(f"{out}.dropmis.imputed.txt",sep="\t", na_rep="NA") 
    return df_imputed

def compute_pca(data):
    if data.shape[1] < 2:  # 只有一个性状时 PC1 贡献率为 100%
        #pc1_df = pd.DataFrame(data.iloc[:, 0], columns=["PC1"])
        pc1_df = pd.DataFrame(data.iloc[:, 0].values, index=data.index, columns=["PC1"])
        return 1.0, pc1_df
    #pca = PCA(n_components=min(data.shape[1], data.shape[0]))  # 维度不能超过样本数
    pca = PCA(n_components=1)  # 只计算第一个主成分
    pc1_scores = pca.fit_transform(data)  # 计算 PC1 得分
    explained_var = pca.explained_variance_ratio_[0]  # PC1 贡献率
    pc1_df = pd.DataFrame(pc1_scores, index=data.index, columns=["PC1"])

    return explained_var,pc1_df


def cor_fig(df,save_path):
    correlation_matrix = df.corr()
    plt.figure(figsize=(10,10))
    sns.clustermap(correlation_matrix,cmap="coolwarm", figsize=(10, 10))
    plt.savefig(f"{save_path}.png", dpi=300)
    correlation_matrix.to_csv(f"{save_path}.cor_matrix.tsv",sep="\t")



if __name__ == "__main__":
    main()
