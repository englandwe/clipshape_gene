#!/usr/bin/env python

#Core classes and functions for the clipshape pipeline

import sys
import os
import random
import re
from Bio import SeqIO,motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from itertools import groupby
from operator import itemgetter

#classes

#handle gtfs with ease and grace
class GtfRec(object):
    def __init__(self,reclist):
         self.seqname=reclist[0]
         self.source=reclist[1]
         self.feature=reclist[2]
         self.start=int(reclist[3])
         self.end=int(reclist[4])
         self.score=reclist[5]
         self.strand=reclist[6]
         self.frame=reclist[7]
         self.attdict={}
         for self.item in reclist[8].strip(';').split('; '):
             self.splitline=self.item.replace('\"','').split(' ')
             self.attdict[self.splitline[0]]=self.splitline[1]

#functions

#understands gencode gtf files
#imports all gtf records with a gene ID (basically everything)
#creates dictionary with gene ids as keys
def importGTFGene(gtffile):
    gtfdict={}
    with open(gtffile) as infile:
        for line in infile:
            if not line.startswith('#'):
                rectmp=GtfRec(line.strip().split('\t'))
                if 'gene_id' in rectmp.attdict.keys():
                    gene=rectmp.attdict['gene_id']
                    try:
                        gtfdict[gene].append(rectmp)
                    except KeyError:
                        gtfdict[gene]=[rectmp]
    return gtfdict

#bit of reformatting to make the generanges look like the txranges file
def outputGeneRanges(gene_id,gtfdict):
    try:
        generec=[rec for rec in gtfdict[gene_id] if rec.feature == 'gene'][0]
        outlist=[gene_id,generec.seqname,generec.strand,generec.start,generec.end]
        return outlist
    except TypeError:
        sys.stderr.write('%s not in gtfdict or lacks gene\n' % gene_id)
        return -1

#exactly what it says
def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final
