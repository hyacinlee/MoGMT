#!/usr/bin/env python3 
import statsmodels.api as sm
import sys 


def main():
    if len(sys.argv) < 3:
        print ("%s  < fpkm file > < Phe file > < Gene bed with 5 Colums >" % (sys.argv[0]))
        exit(1)
    
    ### reads Phe and FPKM data 
    P,Traits,Tsample = readData(sys.argv[2],"h")
    G,Gsample,Genes  = readData(sys.argv[1],"v")

    Pos={}
    ### read Pos 
    for l in open(sys.argv[3],"r"):
        ls=l.strip().split()
        mid = int((int(ls[1])+int(ls[2]))/2)
        Pos[ls[4]] = [ls[0],mid]

    print ("## Total Traits: %s \n\n## Total Genes: %s\n\n## Run Association ... ..." % (len(Traits),len(Genes)))

    cut = 1.0/(len(Genes))
    with open("./cutline.txt","w") as cf:
        cf.write(str(cut))

    for Trait in Traits:
        MIS=0
        ALL=0
        for s in P[Trait]:
            if P[Trait][s] == "NA":
                MIS += 1
            ALL +=1

        if float(MIS)/ALL > 0.25:
            print ("## Trait %s is missing more than 0.25, skipping !!" % (Trait) )
            continue
 
        run_lineRegress(P[Trait],G,Genes,Pos,cut,"Associate.%s" %(Trait))

  



def run_lineRegress(Phe,Ginfo,Genes,Pos,cut,outfile):
    
    of1=open(outfile+".txt","w")
    of1.write("Gene\tChr\tPos\tcoefficient\tcPvalue\tintercept\tiPvalue\n")
    of2=open(outfile+".rdata","w")
    of2.write("SNP\tChromosome\tPosition\tPvalue\n")
    of3=open(outfile+".signal.txt","w")
    of3.write("Gene\tChr\tPos\tcoefficient\tcPvalue\tintercept\tiPvalue\n")
    

    for gene in Genes:
        x=[]
        y=[]
        for sam in Phe.keys():
            if Phe[sam] == "NA" or sam not in Ginfo[gene]:
                continue
            else:
                x.append(float(Ginfo[gene][sam]))
                y.append(float(Phe[sam]))
              
        x = sm.add_constant(x)
        model = sm.OLS(y,x).fit()
        of1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,Pos[gene][0],Pos[gene][1],model.params[1],model.pvalues[1],model.params[0],model.pvalues[0]))
        pp = model.pvalues[1] if not str(model.pvalues[1]) == "nan" else  1.0 
        of2.write("%s\t%s\t%s\t%s\n" % (gene,Pos[gene][0],Pos[gene][1],pp))
        if float(model.pvalues[1]) <= cut:
            of3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,Pos[gene][0],Pos[gene][1],model.params[1],model.pvalues[1],model.params[0],model.pvalues[0]))

        #print (gene,len(x),model.params,model.pvalues,model.params[0],model.params[1],model.pvalues[0],model.pvalues[1])  
        #print (model.summary())
    of1.close()
    of2.close()
    of3.close()

def readData(infile,order):

    P={}
    Hlist=[]
    Vlist=[]

    for l in open(infile,"r"):
        if "t_name" in l:
            continue

        ls = l.strip().split()
        name=ls.pop(0)

        if len(Hlist) ==0:
            Hlist=ls
        else:
            Vlist.append(name)
            for i in range(len(Hlist)):
                k1 = Hlist[i]
                k2 = name
                v = ls[i]
                if order == "h":  ##  H-line  as key
                    P=accumDict(P,k1,k2,v) 
                if order == "v":  ##  V-line  as key
                    P=accumDict(P,k2,k1,v)

    return P,Hlist,Vlist


def accumDict(di,k1,k2,v):
    if k1 not in di:
        di[k1] = {}
    else:
        di[k1][k2] = v 

    return di




main()


