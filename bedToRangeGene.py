#!/usr/bin/env python

#grab an eclip bedfile and extract the relevant bits

import sys
import clipshape_core_gene as csc

def tlbed(bedfile):
    outlist=[]
    with open(bedfile) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            outline=tmpline[0:3]+[tmpline[5]]
            #make 1-based
            outline[1]=str(int(outline[1])+1)
            outline[2]=str(int(outline[2])+1)
            outlist.append(outline)
    return outlist

bedfile=sys.argv[1]
print(csc.flattenList(tlbed(bedfile)))
