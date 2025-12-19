# MoGMT
Multi-omics genetic mapping tools

## Install
 This toolkit is entirely developed based on python. 
 Please ensure that your python >= 3.5 and install the following packages:


Then :
```
pip3 install package.lst 


cd your_path
git clone xxxxxx.git

export PATH="your_path/MoGMT/:$PATH"
export PYTHONPATH="your_path/MoGMT/src:$PYTHONPATH"

MoGMT.py -h
```


## Module Introduction
* **Vcftools**             :   The module for processing vcf file , including format conversion, adding or simplifying information, renaming chromosomes or sample names, etc.
* **Phenotype**         :   The module for processing phenotypic data, including clustering, dimensionality reduction, etc. 
* **Emmax**               :   The module for conducting association analysis using Emmax.
* **Evidence**           :   The module for evidence genes use multi-omics results.
* **Visualize**            :   The module for output Manhantun or other figures.
* **Colocalization**    : The module for perform Co-localization with mutiple traits.
* **BulkSegAna**       :   The module for perform Bulked Segregant Analysis use vcf files.
* **TransAsso**          :   The module for perform Association analysis between gene expression matrix and phenotype in population.


##  Phenotype 

#### 基础操作
```
提取子集(按行或者列)：
MoGMT.py Phenotype -p merge.FPKM.gene.xls --subLines Candidate.co_loc10.ids --out Candidate.co_loc10.FPKM
------------------------------------------------------------

按行统计：
MoGMT.py Phenotype -p Candidate.co_loc10.FPKM.sub_set.tsv  --stat --out Candidate.co_loc10
------------------------------------------------------------

合并combine_files 到 base file
MoGMT.py Phenotype --combine Sample --combine_method Base -p Group1.sub_set.tsv --combine_files Sample.296.SubGroup.txt --out Test.group
------------------------------------------------------------

按列的group统计均值：
MoGMT.py Phenotype -p gene.expression.xls --mean group.txt --out  

head group.txt 
S21_1	FudingDahao
S21_2	FudingDahao
S49_1	Tieluohan
S49_2	Tieluohan
...

------------------------------------------------------------

按列的group进行t-test：
MoGMT.py Phenotype -p gene.expression.grouped.tsv --group_ttest group.txt --out gene.expression

head group.txt
FudingDahao	H
Tieluohan	H
Aijiaowulong L
Echa1	L
...
------------------------------------------------------------


```


#### 聚类
```


```

---
## GWAS 
```
# 基础的 GWAS
MoGMT.py Emmax -v ./snp.dep5.mis90.maf5.simply.vcf -c PCA.cov -p phe.txt

# 如果SNP太多(>2000w)，使用并行提高效率
MoGMT.py Emmax -v ./snp.dep5.mis90.maf5.simply.vcf -c PCA.cov -p phe.txt --split chrom.lst

# 完成后，换阈值及重新筛选
rm -rf GWAS.outputs.done  GWAS.manhantun.done
MoGMT.py Emmax -v ./snp.dep5.mis90.maf5.simply.vcf -c PCA.cov -p phe.txt --split chrom.lst -s 0.000001 

# 示例文件（顺序并不重要）
head PCA.cov 
Sample	PC1	PC2	PC3
10	0.067394	0.0363772	0.023049
104	0.0326593	0.0115734	-0.0513225
105	0.0652511	0.061142	0.0384592
108	0.0467105	0.0194756	-0.017885
110	0.0614426	0.0341291	0.0129854

head chrom.lst
lg1
lg2
lg3
lg4
lg5
```



---
## Evidence
cmd:
```
MoGMT.py Evidence -i ../03.Gwas_singleTriat/Traits.Pos2.emmax.ps.significant.tsv -p ../03.Gwas_singleTriat/Traits.Pos2.phe -o Traits.Pos2.200k -a ../../../../../02.Genotype/02.Genotype_consGenome/03.Annovar/snp.dep5.mis90.maf5.simply -f /work/minghui/Tea_reseq/3.Reseq_FDDB/00.data/refgenome/ref_one/Function.swissport.txt -b /work/minghui/Tea_reseq/3.Reseq_FDDB/00.data/refgenome/ref_one/FD.mRNA.bed -v ../../../4.Area.Grouped.cons/data/snp.dep5.mis90.maf5.simply.vcf.gz --regions 1 200 --regions_score 2 1  -e ../../../../00.data/Evidence/filter.eqtl.cis.txt -w ../../../../00.data/Evidence/weight.txt
```

#### Weight file
```
#Name   Type    Weight  File 
eQTL    eqtl    5       ./add/187_filtered_exp_matrix_qqnorm_file_sort_snp_matrixEQTLoutput_signE5.txt
PopFpkm express 5       ./add/DASZ.220.exp.sorted.FPKM.head.txt
SVs     loci    3       ./add/More2.cluster.svs.bed
trans1  gene    2       ./add/low_vs_high.DEGlist.txt
YL_cmp781       gene    1       ./add/Traits.YL_cmp781.sign.BasicGene.xls
YL_cmp828       gene    1       ./add/Traits.YL_cmp828.sign.BasicGene.xls
```

---
## Visualize

#### plot manhantun
```
MoGMT.py Visualize -m manha --site_value ../03.Gwas_singleTriat/Traits.Pos2.emmax.ps.tsv --hightlight candidate.genes.heightlight --out Traits.Pos2.heightlight


cat heightlight file:
lg12_16125300	PPOX2	red
lg12_35022415	Cyp450	red
lg13_193074629	LAP2	red
lg3_83862148	AAP3	red
```


#### plot boxplot
```
MoGMT.py Visualize -m stat --df Test.group.combine.tsv --box -o Test.group.combine --xData Group --yData Pos2 --ttest CSA1,CSS1 CSA1,CSS2 CSA2,CSS2 CSA2,CSS1 --figsize 8 8 --group_color Group.color

cat Group.color 
CSS1	#3274a1
CSS2	#e1812c
CSA1	#FB5607
CSA2	#8338EC


```


#### plot Cor
```
MoGMT.py Visualize -m stat --df Test.group.combine.tsv --cor -o Test.group.combine --figsize 8 8 --exc Group --plot Row
```



---
### Candidate
```

cat Coloc.txt
Pos2	../../03.Gwas_singleTriat/Traits.Pos2.500k.sign.WeightGenes.xls
Pos5	../../03.Gwas_singleTriat/Traits.Pos5.500k.sign.WeightGenes.xls
Pos15	../../03.Gwas_singleTriat/Traits.Pos15.500k.sign.WeightGenes.xls
Pos35	../../03.Gwas_singleTriat/Traits.Pos35.500k.sign.WeightGenes.xls
Pos53	../../03.Gwas_singleTriat/Traits.Pos53.500k.sign.WeightGenes.xls

MoGMT.py Candidate  -i Coloc.txt --co_local  -o Trait27.express --fit_exp "G128!=None"

```
