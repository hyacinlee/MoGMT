#!/usr/bin/env python3 
import sys
import Common
import argparse
import numpy as np
import pandas as pd
import operator



def get_args(args):
    parser = argparse.ArgumentParser(description="Co localization with mutiple traits.",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-i","--input",required=True,type=str, help="input file, May be TSV or LIST File")
    parser.add_argument("--co_local",action="store_true",default=None, help="Co localization with mutiple traits, then --input must be a list file: Name | File path")
    parser.add_argument("--fit_exp", type=str, help="filter express when reading files,'{HEADER}{Operator}{VALUE}' like 'Express=None'",nargs='+')
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="Candidate") 
    parsed_args = parser.parse_args(args)
    
    return parsed_args
    


def main(args=None):

    args=get_args(args)

    if args.co_local:
        main_co_local(args)  


def main_co_local(args):        

    traits=Common.read_file(args.input,mode="dict",keys=[0],vals=[1])
  
    #print(traits)
    df_list = []
    for t in traits.keys():
        tf = traits[t]
        
        if not Common.check_path_exists(tf,"w"):
            continue

        dt=pd.read_csv(tf,sep="\t")
        if args.fit_exp:
            dt=filter_table(dt,args.fit_exp)
        dt.insert(0, 'Triat', t)
        #print(dt)

        df_list.append(dt)

    df = pd.concat(df_list, ignore_index=True)
    #print(df)

    df.to_csv(f"{args.out}.co_localization.data.xls",sep="\t",index=False)  

    col_start = df.columns.get_loc("Nearest_SNP_P") + 1
    func_cols = df.columns[col_start:]

    gene_info = df.groupby("Gene")[list(func_cols)].first()

    trait_summary = df.groupby("Gene").apply(summarize_traits)

    result = (
        trait_summary
        .join(gene_info)  # 同为 Gene 索引
        .reset_index()
        .sort_values("Total_Base_score", ascending=False)
    )
    #print(result)
    result.to_csv(f"{args.out}.co_localization.merge.xls",sep="\t",index=False)  



def summarize_traits(group):
    triats = sorted(set(group["Triat"]))
    result = {
        "Total_Base_score": group["Base_score"].sum(),
        "Triat_Count": len(triats),
        "Triats": ",".join(triats)
    }

    if "Triat_annot" in group.columns:
        # 构建 Triat -> Triat_annot 映射
        triat_map = dict(zip(group["Triat"], group["Triat_annot"]))
        annot_list = [triat_map[t] for t in triats]
        result["Triat_annots"] = ",".join(annot_list)

    return pd.Series(result)


def filter_table(table: pd.DataFrame, fil_exp: list, mode: str = 'and'):
    """
    过滤 DataFrame。
    fil_exp: list of expressions like ['Express=None', 'Count<20', 'Gene*ABC', 'G128!=None']
    mode: 'and' (默认) 表示所有条件同时满足；'or' 表示任一满足。
    支持操作符: =, !=, <, >, <=, >=, * (包含)
    对于 "None" 值，视为缺失（NaN）或字符串 'None','nan','NA' 或空字符串。
    """
    if mode not in ('and', 'or'):
        raise ValueError("mode 必须是 'and' 或 'or'")

    ops = {
        '=': operator.eq,
        '!=': operator.ne,
        '<': operator.lt,
        '>': operator.gt,
        '<=': operator.le,
        '>=': operator.ge,
        '*': lambda a, b: a.astype(str).str.contains(str(b), na=False)
    }

    df = table.copy()

    # helper: 判断哪些元素应该视为缺失（True = 缺失）
    def missing_mask(series: pd.Series):
        # 将值转为字符串并小写，然后判断是否属于缺失标记集；同时保留 pandas 的 isna 检查
        str_series = series.astype(str).str.strip().str.lower()
        missing_tokens = {'none', 'nan', 'na', ''}
        return series.isna() | str_series.isin(missing_tokens)

    masks = []  # 存放每个表达式对应的 mask

    for exp in fil_exp:
        # 找到操作符（按长度降序，避免把 '!=' 误判成 '!'+'='）
        matched_op = None
        for op_symbol in sorted(ops.keys(), key=len, reverse=True):
            if op_symbol in exp:
                matched_op = op_symbol
                break
        if matched_op is None:
            raise ValueError(f"无效过滤表达式: {exp}")

        header, value = exp.split(matched_op, 1)
        header = header.strip()
        value = value.strip()

        if header not in df.columns:
            raise KeyError(f"列名 '{header}' 不在 DataFrame 中。")

        # 处理 None（用户写 'None'、'none'）这类语义上的缺失
        if value.lower() == 'none':
            if matched_op == '=':
                mask = missing_mask(df[header])
            elif matched_op == '!=':
                mask = ~missing_mask(df[header])
            else:
                # 例如 Count < None 这是不合理的
                raise ValueError(f"不能对 None 使用操作符 {matched_op}: {exp}")
        else:
            # 先尝试把 value 解析为数值（int/float），否则保持为字符串
            parsed_value = value
            is_number = False
            try:
                if '.' in value or 'e' in value.lower():
                    parsed_value = float(value)
                else:
                    parsed_value = int(value)
                is_number = True
            except Exception:
                # 不是整数或浮点数，保留为字符串
                parsed_value = value

            # 包含操作符
            if matched_op == '*':
                mask = ops[matched_op](df[header], parsed_value)
            else:
                # 对于数值比较，尝试把列转换为数值（但不覆盖原列）
                if is_number:
                    # 使用 to_numeric，保留无法转换的为 NaN（在比较时会被处理）
                    col_vals = pd.to_numeric(df[header], errors='coerce')
                    # 当列值为 NaN 且用户比较是数字时，比较结果为 False（自然行为）
                    mask = ops[matched_op](col_vals, parsed_value)
                else:
                    # 字符串/其他比较，直接用 pandas 的比较（astype(str) 保证类型一致）
                    left = df[header].astype(str)
                    # 对 '=' 和 '!=' 做字符串相等比较
                    if matched_op in ('=', '!='):
                        if matched_op == '=':
                            mask = left == str(parsed_value)
                        else:
                            mask = left != str(parsed_value)
                    else:
                        # 对于 <, >, <=, >= 在字符串上比较（按字典序）
                        mask = ops[matched_op](left, str(parsed_value))

        masks.append(mask)

    # 根据 mode 合并 masks
    if not masks:
        return df
    if mode == 'and':
        final_mask = masks[0]
        for m in masks[1:]:
            final_mask = final_mask & m
    else:  # 'or'
        final_mask = masks[0]
        for m in masks[1:]:
            final_mask = final_mask | m

    return df[final_mask]


if __name__ == "__main__":
    main()