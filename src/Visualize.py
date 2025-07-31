#!/usr/bin/env python3 
import sys 
import Common
import argparse
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt


#Chrom   Pos     ID      Value   Order
def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )

    #model_group = parser.add_argument_group('Main optional arguments')

    parser.add_argument("-m","--model",required=True,help="Plot model: manha, pvg, groupGT",choices=["manha", "pvg", "box"],default=None)
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out")
    parser.add_argument("-p","--phe",help="Phenotype file",type=str,default=None)
    parser.add_argument("-f","--figsize",help="Phenotype file",type=str,default=None)
    parser.add_argument("-d","--dpi",help="DPI for png outfiles",type=int,default=300)
    parser.add_argument("-c","--colors",help="colors in fig, one color per line in color file",type=str)

    manha_group = parser.add_argument_group('Manhantun optional arguments (manha)')
    manha_group.add_argument("--site_value",help="Gwas or Other site-value result file")
    manha_group.add_argument("--chrom",help="Chromosom name column",type=str,default="Chrom") 
    manha_group.add_argument("--pos",help="Maker pos column",type=str,default="Pos")
    manha_group.add_argument("--value",help="Value column",type=str,default="Value")
    manha_group.add_argument("--name",help="SNP name column",type=str,default="ID")  
    manha_group.add_argument("--sign",help="Significance cut line: B = 1/n || F = 0.05/n || top1 = top 0.01 || top5 = top 0.05  || a self-defined float",nargs="+",default=["B"])
    manha_group.add_argument("--log",help="change to log10 value",action="store_true",default=True)
    manha_group.add_argument("--sub",help="Draw only one chr",type=str,default=None)
    manha_group.add_argument("--hightlight",help="Hightlight maker ID in files, ids must be same with input",type=str,default=None)
    
    gvp_group = parser.add_argument_group('Correlation between genotype and phenotype optional arguments (pvg)')
    gvp_group.add_argument("--snps",help="Hightlight maker ID in files, ids must be same with input",type=str,default=None)
    gvp_group.add_argument("--gt",help="Traits assosciated GT file( *BasicGT.lst )",type=str,default=None)
    gvp_group.add_argument("--gt_start",help="The start col of GT file",type=int,default=6)
    gvp_group.add_argument("--gt_index",help="The col of maker ID name",type=int,default=1)
    gvp_group.add_argument("--group",help="sample in different group",type=str,default=None)

    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args
 

def main(args=None):

    args=get_args(args)

    if args.model == "pvg":
        PvgMain(args)

    if args.model == "manha":
        ManhantanMain(args)

    #if args.model == "bar":
    #    df = pd.read_table(args.df, sep="\t",)




def ManhantanMain(args):
    
    name_col= {args.chrom:"Chrom",args.pos:"Pos",args.value:"Value",args.name:"ID"}
    df = pd.read_table(args.site_value, sep="\t")
    df = df.rename(columns=name_col)
    df = df.dropna(how="any", axis=0)
  
    df,cut,yname = clean_logs_cut(df,args.log,args.sign,args.value)

    colors=["#3274a1", "#e1812c"] if not args.colors else Common.read_file(args.colors,mode="list",vals=[0],skipAnnot=False)
    print(f"# Manhantun chromosomes color: {colors}")
    hightlight = None if not args.hightlight else Common.read_file(args.hightlight,mode="dict",keys=[0],vals=[1,2])
    print(f"# Hightlight SNP markers: {hightlight} ")

    manhatan_fig_v1(df,cut,yname,args.out,args.sub,hightlight,colors)

 
def manhatan_fig_v1(df, cut, ylab, out, sub=None, hightlight=None,colors=["#3274a1", "#e1812c"],dpi=300):
    """
    绘制曼哈顿图，支持高亮SNP注释、染色体子集和多色，支持双向绘制（根据Order列）。
    df: 包含 Chrom, Pos, Value, ID 四列的数据框（如果有Order列则用于双向绘图）
    cut: 显著性阈值线，可以是一个或多个值的数组（会自动镜像到负方向）
    ylab: y轴标签
    out: 输出文件名（不带扩展名）
    sub: 若指定，只画某条染色体
    hightlight: 字典，key为SNP ID，value为注释文本
    colors: 颜色数组用于交替染色
    dpi: 输出图像分辨率
    """
    df = df.copy()
    
    # 检查是否有Order列，用于双向绘图
    has_direction = 'Order' in df.columns
    
    # 如果指定了子集染色体
    if sub:
        df = df[df['Chrom'] == sub]

    # 排序染色体
    df['Chrom'] = df['Chrom'].astype(str)
    chrom_order = Common.get_sorted_chromosomes(df['Chrom'].unique())
    df['Chrom'] = pd.Categorical(df['Chrom'], categories=chrom_order, ordered=True)
    df = df.sort_values(['Chrom', 'Pos'])

    # 如果有方向信息，调整Value值
    if has_direction:
        df['Value'] = df.apply(lambda row: -row['Value'] if row['Order'] == '-' else row['Value'], axis=1)

    # 计算位置偏移
    df['ind'] = range(len(df))
    df_grouped = df.groupby('Chrom')

    # 绘图
    fig, ax = plt.subplots(figsize=(16, 8))

    x_labels = []
    x_labels_pos = []
    
    for i, (name, group) in enumerate(df_grouped):
        print(f"# Plotting Chrom: {name} into manhattan figure")
        ax.scatter(group['ind'], group['Value'],
                   color=colors[i % len(colors)], s=10)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[0] + group['ind'].iloc[-1]) / 2)

    # 添加cut线（双向）
    for c in cut:
        ax.axhline(y=c, color='gray', linestyle='--', linewidth=1)
        if has_direction:
            ax.axhline(y=-c, color='gray', linestyle='--', linewidth=1)

    # 高亮 SNP 和注释
    if hightlight:
        highlight_ids = list(hightlight.keys())
        highlight_df = df[df['ID'].isin(highlight_ids)]
    
        # 单独处理每个高亮 SNP
        for _, row in highlight_df.iterrows():
            snp_id = row['ID']
            x = row['ind']
            y = row['Value']
            label, color = hightlight.get(snp_id, [snp_id, 'red'])
            
            # 画点（总是改变颜色）
            ax.scatter(x, y, color=color, s=25, zorder=10)
            
            # 只有当注释文本不是"-"时才添加标注
            if label != "-":
                y_text = y + (0.8 if y >= 0 else -0.8)  # 文本高度（根据方向调整）
                ax.annotate(label,
                           xy=(x, y), xycoords='data',
                           xytext=(x, y_text), textcoords='data',
                           arrowprops=dict(arrowstyle='-', color=color, lw=0.8),
                           ha='center', va=('bottom' if y >= 0 else 'top'), 
                           fontsize=8, color=color, rotation=45)

    # 设置轴
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylab)
    ax.set_title(out)
    
    # 如果有方向信息，调整y轴范围并添加零线
    if has_direction:
        y_max = max(df['Value'].max(), abs(df['Value'].min()))
        ax.set_ylim(-y_max*1.1, y_max*1.1)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)


    ticks = ax.get_yticks()
    ax.set_yticklabels([f"{abs(tick):g}" for tick in ticks])

    plt.tight_layout()
    plt.savefig(f"{out}.manhattan.png", dpi=dpi)
    print(f"# Finish output manhattan plot in {out}.manhattan.png")



def manhatan_fig_v0(df, cut, ylab, out, sub=None, hightlight=None,colors=["#3274a1", "#e1812c"],dpi=300):
    """
    绘制曼哈顿图，支持高亮SNP注释、染色体子集和多色。
    df: 包含 Chrom, Pos, Value, ID 四列的数据框
    cut: 显著性阈值线，可以是一个或多个值的数组
    ylab: y轴标签
    out: 输出文件名（不带扩展名）
    sub: 若指定，只画某条染色体
    hightlight: 字典，key为SNP ID，value为注释文本
    colors: 颜色数组用于交替染色
    """
    df = df.copy()
   
    # 如果指定了子集染色体
    if sub:
        df = df[df['Chrom'] == sub]

    # 排序染色体
    df['Chrom'] = df['Chrom'].astype(str)
    chrom_order = sorted(df['Chrom'].unique(), key=chrom_sort_key)
    df['Chrom'] = pd.Categorical(df['Chrom'], categories=chrom_order, ordered=True)
    df = df.sort_values(['Chrom', 'Pos'])

    # 计算位置偏移
    df['ind'] = range(len(df))
    df_grouped = df.groupby('Chrom')

    # 绘图
    fig, ax = plt.subplots(figsize=(16, 8))

    x_labels = []
    x_labels_pos = []
    
    for i, (name, group) in enumerate(df_grouped):
        print(f"# Plotting Chrom: {name} into manhandun figure")
        ax.scatter(group['ind'], group['Value'],
                   color=colors[i % len(colors)], s=10)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[0] + group['ind'].iloc[-1]) / 2)

    # 添加cut线
    for c in cut:
        ax.axhline(y=c, color='gray', linestyle='--', linewidth=1)

    # 高亮 SNP 和注释
    if hightlight:
        highlight_ids = list(hightlight.keys())
        highlight_df = df[df['ID'].isin(highlight_ids)]
    
        # 单独处理每个高亮 SNP（因颜色不同）
        for _, row in highlight_df.iterrows():
            snp_id = row['ID']
            x = row['ind']
            y = row['Value']
            y_text = y + 0.8  # 文本高度
            label, color = hightlight.get(snp_id, [snp_id, 'red'])
    
            # 画点
            ax.scatter(x, y, color=color, s=25, zorder=10)
    
            # 画注释和线
            ax.annotate(label,
                        xy=(x, y), xycoords='data',
                        xytext=(x, y_text), textcoords='data',
                        arrowprops=dict(arrowstyle='-', color=color, lw=0.8),
                        ha='center', va='bottom', fontsize=8, color=color, rotation=45)

    # 设置轴
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylab)
    ax.set_title(out)

    plt.tight_layout()
    plt.savefig(f"{out}.manhantun.png",dpi=dpi)
    print(f"# Finish ouptut manhantun plot in {out}.manhantun.png")


def clean_logs_cut(df,log,sign,name):
    cut = Common.judge_significant(sign,df)  
    
    yname="-log10(P)"
    cut = -np.log10(cut)

    if not log:
        df["P"]=1.0/10**df["P"]
        yname=name
        cut = 1/10**float(cut)
    else:
        df['Value'] = -np.log10(df['Value'])

    return df,cut,yname

def chrom_sort_key(chrom):
    # 提取数字部分作为排序主键；X/Y/MT等排在最后
    try:
        return int(chrom)
    except:
        # 自定义顺序：X=23, Y=24, MT=25，其他未知为99
        chrom_map = {'X': 23, 'Y': 24, 'MT': 25}
        return chrom_map.get(chrom.upper(), 99)



def PvgMain(args):
    
    GTs,header = Common.read_file(args.gt,mode="dict",keys=[args.gt_index-1],vals=[],header=True)
    snps = Common.read_file(args.snps,mode="list",vals=[0])
    Phe = Common.read_file(args.phe,mode="dict",keys=[0],vals=[2])
   
    df_inf = pd.DataFrame.from_dict(Phe, orient='index', columns=['Phenotype'])
    df_inf.index.name="Sample"
    
    if args.group:
        sample_group= Common.read_file(args.group,mode="dict",keys=[0],vals=[1])
        df_inf=Common.updata_df(df_inf,"Group",sample_group,missing="None")
        df_draw=df_inf[df_inf["Group"] != 'None']
        plot_boxplot_ttest(df_draw,'Group', 'Phenotype',"Phe_by_Group.pdf")

    print(df_inf)  

    for snp in snps:
        sample_gt = dict(zip(header[args.gt_start-1:],GTs[snp][args.gt_start-1:]))
        df_inf = Common.updata_df(df_inf,snp,sample_gt,missing="None")
        df_draw1=df_inf[df_inf[snp] != 'NN']
        plot_boxplot_ttest(df_draw1,snp, 'Phenotype',f"{snp}.Phe_by_GT.pdf")
        if args.group:
            df_draw2=df_draw1[df_draw1['Group'] != 'None']
            plot_stacked_bar(df_draw2,f"{snp}.GT_by_group.pdf",value=snp)

    df_inf.to_csv(f"Sample.phe.gt.info.xls",sep="\t",index=True)


def plot_stacked_bar(df,out,value="Phenotype"):
    """
    绘制按照 Group 划分的 Value 堆叠柱状图。
    
    参数:
    df : pd.DataFrame
        包含三列：Sample、Group、Value
    """
    # 确保列名一致
    if not set(['Group', value]).issubset(df.columns):
        raise ValueError(f"DataFrame 必须包含 'Group' 和 {value}三列")

    # 将 Value 转换为计数：每个 Group 内的 Value 出现次数
    count_df = df.groupby(['Group', value]).size().unstack(fill_value=0)

    # 绘图
    ax = count_df.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab20')

    plt.title(f'{value} by Group')
    plt.xlabel('Group')
    plt.ylabel('Count')
    plt.legend(title=value, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(out, format='pdf')


    #print(header)

def plot_boxplot_ttest(df,x,y,out):
    import itertools
    from scipy.stats import ttest_ind

    assert x in df.columns and y in df.columns, "DataFrame must contain 'GT' and 'Phe' columns"
    df = df.dropna(subset=[x, y]).copy()
    df[x] = df[x].astype(str)
    df[y] = pd.to_numeric(df[y], errors='coerce')
    #print(df)
    gt_list = sorted(df[x].unique())
    
    # 箱线图
    plt.figure(figsize=(6, 5))
    ax = sns.boxplot(x=x, y=y, data=df, order=gt_list)
    plt.title(f"{x} vs {y}")

    #if compare==True:
    if len(gt_list) > 1 and len(gt_list) <= 3: # 添加显著性比较（ttest）
        print("# Run t-test compare")
        ymax = df[y].max()
        h = (df[y].max() - df[y].min()) * 0.05
        curr_y = ymax + h
        x_pos = {gt: i for i, gt in enumerate(gt_list)}
        for (gt1, gt2) in itertools.combinations(gt_list, 2):
            group1 = df[df[x] == gt1][y]
            group2 = df[df[x] == gt2][y]
            stat, p = ttest_ind(group1, group2, equal_var=False)
            label = get_sig_label(p)
            x1, x2 = x_pos[gt1], x_pos[gt2]
            ax.plot([x1, x1, x2, x2], [curr_y, curr_y+h, curr_y+h, curr_y], lw=1.5, c='k')
            ax.text((x1 + x2) / 2, curr_y + h + 0.01, label, ha='center', va='bottom')
            curr_y += h * 2

    plt.tight_layout()

    plt.savefig(out, format='pdf')
    plt.close()


def get_sig_label(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'


if __name__ == '__main__':

    main()
