#!/usr/bin/env python3 
import sys
import Common
import random
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.stats import ttest_ind
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
    parser.add_argument("-b","--bootstrap",help="Bootstrap times to run Cluster",default=500,type=int)
    parser.add_argument("--index",help="Add or change Index Name in output file",default="Sample",type=str)
    parser.add_argument("--mean",help="Calculate the average value of multiple columns using the specified grouping： name | groups",default=None,type=str)
    parser.add_argument("--group_ttest",help="Conduct a t-test using the specified group: name | groups",default=None,type=str)
    parser.add_argument("--stat",help="",default=None,action='store_true')
    parser.add_argument("--combine",help="Combine Column names",default=None)
    parser.add_argument("--combine_files",help="Files combined to -phe",default=None,nargs="*")
    parser.add_argument("--combine_method",help="Files combined method,choose from Base, Union or Intersection",default="Base")
    #parser.add_argument("--rename",help="Replace the column names with the specified name relationship：  old name | new name",default=None,type=str)
    parser.add_argument("--subColumns",help="Extract the Columns in the file, in columns.",default=None,type=str)
    parser.add_argument("--subLines",help="Extract the Lines in the file, in lines.",default=None,type=str)
    parser.add_argument("--transpose",help="transpose line and columns and output",action='store_true')
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out")   
    parsed_args = parser.parse_args(args)
    
    return parsed_args
    

#corr_threshold=0.8, pc_threshold=0.8

def main(args=None):

    args=get_args(args)

    phe_data=load_data(args.phe)

    if args.mean:
        run_matrix_mean(phe_data,args)

    if args.group_ttest:
        run_matrix_group_ttest(phe_data,args)

    elif args.subColumns or args.subLines:
        run_matrix_sub(phe_data,args)

    elif args.transpose:
        run_matrix_transpose(phe_data,args)

    elif args.stat:
        row_statistics(phe_data,args.out)

    elif args.combine:
        lis=[args.phe]+args.combine_files
        run_combine(file_list=lis, on_col=args.combine, output=f"{args.out}.combine.tsv", method=args.combine_method)
    else:
        run_phe(phe_data,args)


def run_combine(file_list,on_col,output,method="Union",fillna="None", sep="\t",):

    if method not in ["Base", "Union", "Intersection"]:
        raise ValueError("method must be one of: Base, Union, Intersection")

    merged = None
    key_sets = []

    # === 第一步：读取所有文件并记录 on_col 值 ===
    dfs = []
    for f in file_list:
        df = pd.read_csv(f, sep=sep, dtype=str)
        if on_col not in df.columns:
            raise ValueError(f"{f} 中缺少用于合并的列: {on_col}")

        dfs.append((f, df))
        key_sets.append(set(df[on_col].astype(str)))

    # === 第二步：计算保留的键 ===
    if method == "Base":
        keep_keys = key_sets[0]  # 只保留第一个文件的 key
    elif method == "Union":
        keep_keys = set.union(*key_sets)
    elif method == "Intersection":
        keep_keys = set.intersection(*key_sets)
    else:
        raise ValueError("Invalid merge method")

    # === 第三步：过滤 df ===
    filtered_dfs = []
    for f, df in dfs:
        df = df[df[on_col].astype(str).isin(keep_keys)].copy()
        filtered_dfs.append(df)

    # === 第四步：依次合并 ===
    merged = filtered_dfs[0]
    for df in filtered_dfs[1:]:
        merged = pd.merge(merged, df, on=on_col, how="outer")

    # === 填充缺失值 ===
    merged = merged.fillna(fillna)

    # === 输出 ===
    if output:
        merged.to_csv(output, sep=sep, index=False)
    else:
        print(merged)

    return merged



def row_statistics(df,output):
    df_num = df.select_dtypes(include=[np.number])

    stats_df = pd.DataFrame({
        'sum': df_num.sum(axis=1),
        'mean': df_num.mean(axis=1),
        'std': df_num.std(axis=1),
        'count': df_num.count(axis=1),
        'nonzero': (df_num != 0).sum(axis=1)
    })
    stats_df = stats_df.reset_index()
    stats_df.rename(columns={'index': 'Info'}, inplace=True)
    stats_df.to_csv(f"{output}.stat_row.tsv",sep="\t",index=False)


def run_matrix_group_ttest(df,args):
    print("Conduct a t-test using the specified group")
    sample_to_group=Common.read_file(args.group_ttest,mode="dict",vals=[1],keys=[0])
    valid_samples = [s for s in df.columns if s in sample_to_group]
    df = df[valid_samples]

    groups = list(set(sample_to_group[s] for s in valid_samples))
    if len(groups) != 2:
        raise ValueError(f"Exactly two sets are needed, but there are now {len(groups)} groups: {groups}")

    g1, g2 = groups
    g1_samples = [s for s in valid_samples if sample_to_group[s] == g1]
    g2_samples = [s for s in valid_samples if sample_to_group[s] == g2]

    results = []
    for gene, row in df.iterrows():
        vals1 = row[g1_samples].astype(float).values
        vals2 = row[g2_samples].astype(float).values
        
        # T-test
        stat, pval = ttest_ind(vals1, vals2, equal_var=False, nan_policy="omit")
        
        g1_mean, g1_std = vals1.mean(), vals1.std(ddof=1)
        g2_mean, g2_std = vals2.mean(), vals2.std(ddof=1)
        
        # log2 fold change
        #log2fc = np.log2((g2_mean + 1) / (g1_mean + 1))
        
        results.append({
            "gene_id": gene,
            "pvalue": pval,
            #"log2FC": log2fc,
            f"{g1}_mean": g1_mean,
            f"{g1}_std": g1_std,
            f"{g2}_mean": g2_mean,
            f"{g2}_std": g2_std     
        })

    res_df = pd.DataFrame(results).sort_values("pvalue")
    res_df.to_csv(f"{args.out}.grouped.t-test.tsv",sep="\t",index=False)
    print(res_df)
    #return res_df



def run_matrix_mean(df,args):
    print("Calculate the average value of multiple columns using the specified grouping")
    sample_to_group=Common.read_file(args.mean,mode="dict",vals=[1],keys=[0])
    valid_samples = [s for s in df.columns if s in sample_to_group]
    df = df[valid_samples]
    df_renamed = df.rename(columns=sample_to_group)
    df_grouped = df_renamed.groupby(level=0, axis=1).mean()
    df_grouped.to_csv(f"{args.out}.grouped.tsv",sep="\t")
    #print(df_grouped)


def run_matrix_sub(phe_data,args):
    print("Get Sub df with samples or traits")
    df_sub=sub_sample_traits(phe_data,args.subColumns,args.subLines)
    df_sub.to_csv(f"{args.out}.sub_set.tsv",sep="\t")
    exit(0)

def run_matrix_transpose(phe_data,args):
    df = phe_data.T
    df.index.name = args.index
    df.to_csv(f"{args.out}.transpose.tsv",sep="\t")
    exit(0)



def run_phe(phe_data,args):

    df_mis = missing_values(phe_data, args.miss_threshold,args.out)

    df_mis_inpute = df_mis
    if args.impute:
        df_mis_inpute=impute_values(df_mis,args.n_neighbors,args.out)


    if args.globals:

        best_group_length = float('inf')
        best_n50 = -1
        best_uncluster = float('inf')
        best_top20 = -1
        result=[]
        change=[]

        for i in range(args.bootstrap):
            groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold)

            group_length=len(groups)

            n50,top_20_ratio,uncluster = n50_from_2d_array(groups)

            is_better = False
    
            if group_length < best_group_length:
                is_better = True
            elif group_length == best_group_length:
                if n50 > best_n50:
                    is_better = True
                elif n50 == best_n50:
                    if uncluster < best_uncluster:
                        is_better = True
                    elif uncluster == best_uncluster:
                        if top_20_ratio > best_top20:
                            is_better = True

            if is_better:
                best_group_length = group_length
                best_n50 = n50
                best_uncluster = uncluster
                best_top20 = top_20_ratio
                result=[groups,pc_values,group_dfs,group_pc_dfs]

            change.append([i+1,best_group_length,best_n50,best_uncluster,best_top20])
            print(f"# Bootstrap {i+1} times，group lengths = {len(groups)} ;N50 = {n50} ; top_20_ratio = {top_20_ratio} ;  uncluster = {uncluster}")

        groups,pc_values,group_dfs,group_pc_dfs=result
        output(df_mis_inpute,groups,pc_values,group_dfs,group_pc_dfs,args.out,args.fig)
        df_change = pd.DataFrame(change, columns=['times', 'length', 'n50', 'uncluster', 'top20'])
        df_change.to_csv(f"{args.out}.best_history.tsv", sep='\t', index=False)


    if args.locals: #locals
        start_list=pd.read_csv(args.locals, sep='\t').iloc[:, 0].tolist()
        groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold,start_list)
        output(df_mis_inpute,groups,pc_values,group_dfs,group_pc_dfs,args.out,args.fig)



    #print(df_mis_inpute)
'''    
    start_list = df_mis_inpute.columns.tolist()
    if args.locals:
        start_list=pd.read_csv(args.locals, sep='\t').iloc[:, 0].tolist()

    shuffled_lists = [random.sample(start_list, len(start_list)) for _ in range(args.bootstrap)]
    #print(shuffled_lists)

    min_len = float('inf')          # 当前最小的 len(groups)
    best_result = None              # 保存对应的最优结果

    for i,order_list in enumerate(shuffled_lists):
        groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold)
        current_len = len(groups)

    # 如果更小，就更新最优结果
        if current_len < min_len:
            min_len = current_len
            best_result = (groups, pc_values, group_dfs, group_pc_dfs, order_list)

    # 每运行5次打印一次进度
        if (i + 1) % 1 == 0:
            print(f"# Bootstrap {i + 1} times，group lengths = {min_len}")
'''

def n50_from_2d_array(arr2d):
    lengths = [len(row) for row in arr2d]
    
    lengths_sorted = sorted(lengths, reverse=True)
    #print(lengths_sorted)
    top_20_percent_count = max(1, int(len(lengths_sorted) * 0.2))  # 至少 1 个
    top_20_sum = sum(lengths_sorted[:top_20_percent_count])
    total_sum = sum(lengths_sorted)
    top_20_ratio = top_20_sum / total_sum

# 2️⃣ 统计长度为 1 的元素个数
    count_len1 = sum(1 for l in lengths if l == 1)

    total = sum(lengths_sorted)
    cum_sum = 0
    
    for l in lengths_sorted:
        cum_sum += l
        if cum_sum >= total / 2:
            n50 = l
            break
    
    return [n50,top_20_ratio,count_len1]

def sub_sample_traits(df,subColumns,subLines):
    #print(df.columns)
    df.index = df.index.astype(str)
    if subColumns:
        trait_list=Common.read_file(subColumns,mode="list",vals=[0])
        print(trait_list)
        missing_traits = [t for t in trait_list if t not in df.columns]
        sub_traits = [t for t in trait_list if t in df.columns]
        if missing_traits:
            print(f"The following traits were not found in the DataFrame and will be ignored: {missing_traits}")
        df = df.loc[:,sub_traits]
    if subLines:
        samples_list=Common.read_file(subLines,mode="list",vals=[0])
        print(samples_list)
        df.index = df.index.astype(str)
        #print(df.index)
        valid_samples = [str(s) for s in samples_list if s in df.index]
        #print(samples_list)
        #print(df)
        print(valid_samples)
        df = df.loc[valid_samples]
        #print(df)
    return df


# 进行性状分组
def group_traits(df, corr_threshold=0.8, pc_threshold=0.8, start_list=None):
    corr_matrix = df.corr()  # 保留正负号

    #all_traits = set(df.columns)
    all_traits = list(df.columns)
    random.shuffle(all_traits)
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
                #print (f"#Cluster starts with {base_trait}")
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
    df = pd.read_csv(file_path,sep="\t", index_col=0,dtype={0: str})  # 以第一列（样本名）作为索引
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
