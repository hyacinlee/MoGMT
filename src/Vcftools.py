#!/usr/bin/env python3
import re
import sys
import Common
import argparse


def get_args(args):
    parser = argparse.ArgumentParser(description='''
    The module for processing vcf file , including format conversion, adding or simplifying information, renaming chromosomes or sample names, etc

    This is merely a functional supplement to vcftools. If you need to perform more professional operations on vcf, please use the original vcftools",
    ''',
    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-v", "--vcf",required=True, help="Path of the input vcf file")
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="Output.vcf")
    parser.add_argument("-c","--chr",help="Keep SNPs of Chromosomes in files, if --rename_chr is activated, two columns will be required",type=str,default=None)
    parser.add_argument("-k","--keep",help="Keep sampels in files, if --rename_samples is activated, two columns will be required",type=str,default=None)
    parser.add_argument("-t","--tags",help="Retain the information of the specified tags like [ GT AD DP ...]",nargs="+",default=None)
    parser.add_argument("-a","--add_ids",help="Add SNP-ID information in the format of chr_pos",action="store_true",default=None)
    parser.add_argument("-d","--diploid",help="Convert the haploid GT to the pseudo diploid GT, 0 -> 0/0, 1 -> 1/1",action="store_true",default=None)
    parser.add_argument("-s","--simplify",help="Filter out all unnecessary information and retain only the simplest genotype.",action="store_true",default=None)
    parser.add_argument("--keep_header",help="Keep vcf header info in simplify mode",action="store_true",default=None)
    parser.add_argument("--rename_chr",help="Rename the Chromosomes according to the relationship specified in --chr. Chromosomes not in the file will be discarded",action="store_true",default=None)
    parser.add_argument("--rename_sample",help="Rename the Sampels according to the relationship specified in --keep. Samples not in the file will be discarded",action="store_true",default=None)


     
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args
  
def main(args=None):

    args=get_args(args)

    keep_chrom = [x[0] for x in Common.read_file(args.chr,mode="list")] if args.chr else []    # keep chroms 
    if args.rename_chr:
        change_chrom = Common.read_file(args.chr, mode="dict", keys=[0] ,vals=[1])
    #print(change_chrom)

    keep_sample = [x[0] for x in Common.read_file(args.keep,mode="list")] if args.keep else []
    #change_sample = dict(zip(keep_sample))
    if args.rename_sample:
        change_sample= Common.read_file(args.keep, mode="dict", keys=[0] ,vals=[1])


    with open(args.vcf) as fin,open(args.out, 'w') as fout:
        for line in fin:
            if line.startswith('##'):      # annovar info
                if args.keep_header:
                    fout.write(line)
                elif args.simplify:
                    continue
                else:
                    fout.write(line)

            elif line.startswith('#CHROM'): # sample info 
                fout.write(line)

            else:
                parts = line.strip().split('\t')

                if args.chr:
                    if parts[0] not in keep_chrom:
                        continue
                    if args.rename_chr:
                        parts[0] = change_chrom[parts[0]]

                if args.add_ids:
                    parts[2] = f"{parts[0]}_{parts[1]}"


                all_tags  = parts[8].split(':')   # all tags 
                need_tags = all_tags

                if args.simplify:
                    parts[7] = "."  
                    parts[6] = "."
                    parts[5] = "."
                    need_tags = ["GT"]

                if args.tags:
                    need_tags = args.tags  

                parts[8] = ":".join(need_tags)


                if not need_tags == all_tags or args.diploid:
                    for i in range(9, len(parts)):
                        parts[i] = format_sample_info(parts[i],all_tags,need_tags,args.diploid)
                
                fout.write("\t".join(parts)+"\n")



def format_sample_info(inf,formats,need_tags,diploid=None):
    #print(diploid) 
    info_dict=dict(zip(formats,inf.split(':')))
    if diploid:
        info_dict["GT"] = info_dict["GT"]+"/"+info_dict["GT"]
    new_inf=(":".join([info_dict[x] for x in need_tags]))
    return new_inf


    #simplify_and_filter_vcf(sys.argv[1], sys.argv[2])
def rename_header():
    pass



def simplify_and_filter_vcf(input_file, output_file):
    with open(input_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            parts = line.strip().split('\t')
            format_fields = parts[8].split(':')

            try:
                gt_index = format_fields.index("GT")
                ad_index = format_fields.index("AD")
                dp_index = format_fields.index("DP")
            except ValueError:
                continue

            new_sample_fields = []
            all_genotypes = set()
            skip_line = False

            for i in range(9, len(parts)):
                fields = parts[i].split(':')
                if len(fields) <= max(gt_index, ad_index, dp_index):
                    new_sample_fields.append('./.:.:.')
                    all_genotypes.add('./.')
                    continue

                raw_gt = fields[gt_index]
                if '.' in raw_gt:
                    new_sample_fields.append('./.:.:.')
                    all_genotypes.add('./.')
                    continue

                alleles = re.split(r'[|/]', raw_gt)
                if len(alleles) != 2:
                    new_sample_fields.append('./.:.:.')
                    all_genotypes.add('./.')
                    continue

                alleles = sorted(alleles)
                new_gt = '/'.join(alleles)

                ad = fields[ad_index]
                dp = fields[dp_index]

                new_sample_fields.append(f"{new_gt}:{ad}:{dp}")
                all_genotypes.add(new_gt)

            if len(all_genotypes) <= 1:
                continue

            parts[8] = "GT:AD:DP"
            parts[9:] = new_sample_fields
            fout.write('\t'.join(parts) + '\n')



if __name__ == "__main__":
    main() 
