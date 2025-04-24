import sys
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from itertools import combinations

# 读取数据
# 1. 读取数据
def load_data(file_path):
    df = pd.read_csv(file_path,sep="\t", index_col=0)  # 以第一列（样本名）作为索引
    df.replace([np.inf, -np.inf], np.nan, inplace=True)  # 处理 inf
    return df

# 2. 处理缺失值（20% 阈值 + KNN 填充）
def handle_missing_values(df, missing_threshold=0.4, n_neighbors=5):
    print(f"数据集形状: {df.shape}")  # 输出数据的行列数
    df = df.select_dtypes(include=[np.number])  # 只保留数值型列
    
    # 计算每列缺失率
    missing_rates = df.isna().mean()
    
    # 删除缺失率超过 `missing_threshold`（默认为 20%）的列
    columns_to_drop = missing_rates[missing_rates > missing_threshold].index
    df.drop(columns=columns_to_drop, inplace=True)
    
    print(f"删除缺失率超过 {missing_threshold*100:.0f}% 的性状列: {list(columns_to_drop)}")

    # 用 KNN 填充剩余缺失值
    imputer = KNNImputer(n_neighbors=n_neighbors)
    df_imputed = pd.DataFrame(imputer.fit_transform(df), columns=df.columns, index=df.index)
    
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
    correlation_matrix = df.T.corr()
    plt.figure(figsize=(10,10))
    sns.clustermap(correlation_matrix,cmap="coolwarm", figsize=(10, 8))
    plt.savefig(save_path, dpi=300)


# 进行性状分组
def group_traits(df, corr_threshold=0.6, pc_threshold=0.8):
    corr_matrix = df.corr().abs()  # 计算相关性矩阵（取绝对值）
    remaining_traits = set(df.columns)  # 未分组的性状
    grouped_traits = []  # 存储最终分组
    group_pc_values = []  # 存储每个组的 PC1 贡献率
    group_dfs = []  # 存储每个组的数据
    group_pc_dfs=[] # 存储每个组的PC1的DF

    while remaining_traits:
        group = []  # 当前新的一组
        base_trait = remaining_traits.pop()  # 选一个基准性状
        group.append(base_trait)

        # 找出所有与 base_trait 相关性 > 0.8 的性状
        candidates = {trait for trait in remaining_traits if corr_matrix.loc[base_trait, trait] > corr_threshold}
        
        temp_group = group.copy()
        for trait in candidates:
            temp_group.append(trait)
            temp_df = df[temp_group].dropna()
            if temp_df.shape[0] > 1 and compute_pca(temp_df)[0] >= pc_threshold:
                group.append(trait)
                remaining_traits.remove(trait)  # 从未分组集合中移除
            else:
                temp_group.pop()  # 如果 PC1 < 80%，移除该性状


        if len(group) == 0:
            continue  # 避免添加空 group

        group_df = df[group].dropna()
        if group_df.shape[1] == 0:
            continue  # 避免空 group 进入最终结果


        # 获取该组最终 PCA 贡献率
        group_df = df[group].dropna()
        group_pc, pc1_scores = compute_pca(group_df) #if group_df.shape[0] > 1 else (1.0, pd.DataFrame(group_df.iloc[:, 0], columns=["PC1"]))
        grouped_traits.append(group)
        group_pc_values.append(group_pc)
        group_pc_dfs.append(pc1_scores)
        group_dfs.append(group_df)


    sorted_groups = sorted(zip(grouped_traits, group_pc_values, group_dfs,group_pc_dfs), key=lambda x: len(x[0]), reverse=True)

    return zip(*sorted_groups)  # 返回拆解后的 (groups, pc_values, group_dfs)

    #return grouped_traits, group_pc_values, group_dfs

# 主函数
def main():
    file_path = sys.argv[1]  # 你的输入文件
    out_path  = "Combine.Triat.txt"
    out=open(out_path,"w")

    df2 = load_data(file_path)
    df2 = handle_missing_values(df2)
    
    groups,pc_values,group_dfs,group_pc_dfs = group_traits(df2)

    df1 = pd.DataFrame(index=df2.index)
    for i, pc1_scores in enumerate(group_pc_dfs, start=1):
        df1[f"Group{i}"] = pc1_scores["PC1"]    

    
    # 输出分组结果
    for i, (group, pc1,dfg,dfpc) in enumerate(zip(groups, pc_values,group_dfs,group_pc_dfs), start=1):
        if len(group) > 2 :
            cor_fig(dfg,"Combine.Group%s.png"%(i))
            #out_path  = "Combine.Triat.txt"                 
        
        #dfpc.to_csv(f"Combine.Group{i}s.pc1",sep="\t")  # PC out      
        print(f"Group{i}: {group}, PC1 贡献率: {pc1:.2%}")        
        out.write(f"Group{i}\t{pc1:.2%}\t{group}\n")     # mian out 


    selected_traits = [trait for group in groups if len(group) >= 2 for trait in group]
    df3 = df2[selected_traits]

    selected_traits_df4 = [group[0] for group in groups if len(group) == 1]
    df4 = df2[selected_traits_df4]

    df1.to_csv("Sample.PC1_scores.txt",sep="\t")  # 每个样本对应每个group的PC1值
    df2.to_csv("Sample.Trait.txt",sep="\t")  # 过滤和填充后的性状矩阵
    df3.to_csv("Sample.Trait.inGroup.txt",sep="\t")  # 仅包含 group 长度 >= 2 的性状
    df4.to_csv("Sample.Trait.unGroup.txt",sep="\t")


if __name__ == "__main__":
    main()

