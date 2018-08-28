#!/usr/bin/env python

#ENST00000425361        2       +       130190577       0.704   NEWVAL 1       0       T       1

import sys

subline=[]
with open(sys.argv[1]) as infile:
    prevline=next(infile).strip().split('\t')
    subline.append(prevline)
    for line in infile:
        tmpline=line.strip().split('\t')
        if tmpline[-1] == prevline[-1]:
            #consecutive
            subline.append(tmpline)
        else:
            #wrap up old range
            if subline[0][2] == '-':
                subline=subline[::-1]
            min_idx=-(len(subline)/2) #yes, floor
            for entry in subline:
                entry.append(min_idx)
                print('\t'.join([str(x) for x in entry]))
                min_idx += 1
            #new range
            subline=[]
            prevline=list(tmpline)
            subline.append(prevline)
#finish last range
    min_idx=-(len(subline)/2) #yes, floor
    for entry in subline:
        entry.append(min_idx)
        print('\t'.join([str(x) for x in entry]))
        min_idx += 1

