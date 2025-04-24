#!/usr/bin/env python3 
import argparse
import importlib
import os
import sys

#from src import Emmax
# python module need
#   numpy
#   pandas
#   seaborn
#   matplotlib
#   sklearn
#   subprocess

SUB_COMMANDS = {
    "Emmax": "Emmax", 
    "Evidence": "Evidence",
    "Phenotype": "Phenotype",
    "Visualize": "Visualize",
    "Vcftools" : "Vcftools",
    "BulkSegAna":"BulkSegAna"
}

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='''
    Main Prgram of Multi-omics genetic mapping tools , Plase choose Function from:

        Vcftools     :   The module for processing vcf file , including format conversion, adding or simplifying information, renaming chromosomes or sample names, etc.

        Phenotype    :   The module for processing phenotypic data, including clustering, dimensionality reduction, etc. 

        Emmax        :   The module for conducting association analysis using Emmax.

        Evidence     :   The module for evidence genes use multi-omics results.

        Visualize    :   The module for output Manhantun or other figures. 

        BulkSegAna   :   The module for perform Bulked Segregant Analysis use vcf files.
        
    Date      : 2025-04-14
    Version   : 0.1
        ''',
        add_help=False,
        usage="%(prog)s <command> [args]"
    )
    parser.add_argument(
        "command",
        choices=SUB_COMMANDS.keys(),
        help="Avaliable Functions",
    )
    
    if "-h" in sys.argv or "--help" in sys.argv:
        if len(sys.argv) == 2:
            parser.print_help()
            sys.exit(0)
        else:
            cmd = sys.argv[1]
            main_func = load_module(cmd)
            main_func(["-h"])  
            sys.exit(0)

    # 正常解析命令名
    args, remaining_args = parser.parse_known_args()
    
    # 加载并执行子命令
    main_func = load_module(args.command)
    main_func(remaining_args)  # 传递剩余参数

def load_module(module_name):
    """动态加载子命令模块并返回其 main 函数"""
    try:
        module = importlib.import_module(SUB_COMMANDS[module_name])
        #print(os.path.exists("src/Emmax.py"))
        #module = importlib.import_module("src.Emmax")
        return module.main
    except ImportError:
        print(f"Error: module {module_name} not exist")
        sys.exit(1)
    except KeyError:
        print(f"Error: unkonw sub-Function '{module_name}'")
        sys.exit(1)


if __name__ == "__main__":
    main()
