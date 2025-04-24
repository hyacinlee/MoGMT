#!/usr/bin/env python 
import sys 

if not  len(sys.argv) == 4:
    print "%s <vcf>  <trait> < *signals.csv >" % (sys.argv[0])
    exit(1)

sites=[]
for line in open(sys.argv[3]):
    if line.startswith("\"SNP\""):
        continue
    #s = line.strip().split(",")[0].replace("\"","")
    s = line.replace("\"","").split(",") 
    #ss == s.split("-")    
    sites.append("%s-%s"%(s[1],s[2]))
#print sites

gt={}
sample=[]
for line in open(sys.argv[1],"r"):
    if line.startswith("##"):
        continue
    s = line.strip().split()
    tag="%s-%s" %(s[0],s[1])
    if line.startswith("#CHROM"):
        sample=s
    elif tag in sites:
        gt[tag] = {}
        for i in range(9,len(s)):
            gg=s[i].split(":")[0]
            gt[tag][sample[i]] = gg
    else:
        if len(gt) == len(sites):
            break
        else:
            continue
#print gt
print "#Pos\tSample\tGT\tTrait"
for site in sites:
    fo=open(sys.argv[2],"r")
    for l in fo:
        s = l.strip().split("\t")
        if l.startswith("Taxa"):
            continue
        else:
            if s[0] in gt[site]:
                print "%s\t%s\t%s\t%s" %(site,s[0],gt[site][s[0]],s[1])
    fo.close()
