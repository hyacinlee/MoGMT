#!/usr/bin/env python 
import sys 

sample=[]
stat={}
for line in open(sys.argv[1]):
    s=line.strip().split()
    if line.startswith("##"):
        continue
    elif line.startswith("#CHR"):
        sample=s[9:]
        #print sample
    else:
        tvts = "transversions"
        if (s[3] == "A" and s[4] == "G") or ( s[3] == "G" and s[4] =="A" ) or (s[3]=="C" and s[4] =="T") or (s[3]=="T" and s[4]=="C"):
            tvts = "transitions"

        gts=s[9:]
        for i in range(len(gts)):
            ss = sample[i]
            if ss not in stat:
                stat[ss]=[0,0,0,0,0,0]   # mis, var ,hom, het, ts, tv
            gt = gts[i].split(":")[0]
            if gt == "./." or gt == ".":
                stat[ss][0] += 1
#            elif gt == "0/0" or gt == "0|0":
#                continue
            else:
                stat[ss][1] += 1
                if tvts == "transitions":
                    stat[ss][4]+=1
                else:
                    stat[ss][5]+=1
   
                if gt == "1/1" or gt == "1|1" or gt == "0/0" or gt == "0|0":
                    stat[ss][2]+=1
                else:
                    stat[ss][3]+=1

print "Sample\tMis\tNumber\tHom\tHet\tTs\tTv"
for sa in sample:
    print "%s\t%s" %(sa,"\t".join(map(str,stat[sa])))



