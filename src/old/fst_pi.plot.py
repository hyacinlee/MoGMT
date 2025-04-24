#!/usr/bin/env python 

import sys 
import pandas as pd
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

def main():
    #infile = "C:/Users/diggers/Desktop/combine_fst_ROD.test.xls"
    infile = sys.argv[1] 
    data = pd.read_table(infile,sep='\t',header=0)
    #data = data.dropna()
    #print data
    bins=100
    #set axe
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    nullfmt = NullFormatter() 
    data.max() 
    
    plt.rcParams['figure.figsize'] = (10,10)
    plt.rcParams['savefig.dpi'] = 300

    ax1 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_scatter)
    ax3 = plt.axes(rect_histy)
    ax1.xaxis.set_major_formatter(nullfmt)
    ax3.yaxis.set_major_formatter(nullfmt)
    
    # data format
    cut_index =  int(0.05*len(data['ROD']))
    cut1l     =  sorted(data['ROD'],reverse=False)[cut_index]  #data['ROD'].rank(ascending = True)[cut_index]
    cut1h     =  sorted(data['ROD'],reverse=True)[cut_index] # data['ROD'].rank(ascending = False)[cut_index]
    cut2      =  sorted(data['FST'],reverse=True)[cut_index]  #data['FST'].rank(ascending = False)[cut_index]
    print "Fst cut-line %s " %(cut2)
    print "Rod cut-line %s and %s" %(cut1l,cut1h)
    #print (data['ROD'].rank(ascending = False))
    d1h = data[ data.ROD >= cut1h ]
    d1m = data[ (data.ROD > cut1l) & (data.ROD < cut1h) ]
    d1l = data[ data.ROD <= cut1l ]
    d2l = data[ data.FST <= cut2 ]
    d2h = data[ data.FST > cut2 ] 
    c3  = colored(data['ROD'],data['FST'],cut1l,cut1h,cut2)
    maxROD = data.max()['ROD']
    maxFST = data.max()['FST']
    print "ROD max is %s " %(maxROD)
    print "FST max is %s " %(maxFST)
    ax1.set(xlim=[0, data.max()['ROD']],ylabel="Frequency(%)")
    ax2.set(xlim=[0, data.max()['ROD']], ylim=[0,data.max()['FST']] ,ylabel='Fst', xlabel='Rod')
    ax3.set(ylim=[0, data.max()['FST']],xlabel="Frequency(%)")

    ax1.hist([d1l['ROD'],d1m['ROD'],d1h['ROD']],bins=bins,stacked=True,density=True,color=["blue","grey",'blue'])
    ax12 = ax1.twinx()
    ax12.set_ylabel('Cumulative(%)')
    ax12.hist(data.ROD, bins=bins, density=True, histtype='step', cumulative=True,color="black")

    ax2.scatter(data['ROD'],data['FST'],c=c3,s=1)

    ax3.hist([d2l['FST'],d2h['FST']],bins=bins,color=["grey","red"],stacked=True,density=True,orientation='horizontal')  #,orientation='horizontal'
    ax32 = ax3.twiny()
    ax32.set_xlabel('Cumulative(%)')
    ax32.hist(data.FST, bins=bins, density=True, histtype='step', cumulative=True,orientation='horizontal',color="black")

    ax1.axvline(x=cut1h,ls="--",c="grey",lw=1)
    ax1.axvline(x=cut1l,ls="--",c="grey",lw=1)
    ax2.axvline(x=cut1h,ls="--",c="grey",lw=1)
    ax2.axvline(x=cut1l,ls="--",c="grey",lw=1)
    ax2.axhline(y=cut2,ls="--",c="grey",lw=1)
    ax3.axhline(y=cut2,ls="--",c="grey",lw=1)
    ax3.axhline(y=cut2,ls="--",c="grey",lw=1)
    #plt.show()
    plt.savefig('%s.png' %(sys.argv[1]), dpi=300)

    out1 = data[ (data.ROD > cut1h) & (data.FST >cut2) ].sort_values(by=["CHR","POS"])
    out2 = data[ (data.ROD < cut1l) & (data.FST >cut2) ].sort_values(by=["CHR","POS"])

    out1.to_csv("%s.hight.candiate" % sys.argv[1], sep ="\t",index=False,header=False)
    out2.to_csv("%s.low.candiate" % sys.argv[1], sep ="\t",index=False,header=False)

def colored(d1,d2,cut1l,cut1h,cut2):
    c3=[]
    for i in range(len(d1)):
        if d2[i] >= cut2:  # fst 
            if d1[i] < cut1l or d1[i] > cut1h:  # rod
                c3.append("green")
            else:
                c3.append("black")
        else:
            c3.append("black")

    return (c3)


if __name__ == "__main__":
    main()
