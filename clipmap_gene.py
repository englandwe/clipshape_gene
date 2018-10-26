#!/usr/bin/env python

#USAGE: shapedata clipdata gtf output_prefix

#returns:  a file of gene ranges (<output_prefix>.generanges)
#returns:  a file of mapped clip & shape data (<output_prefix>.clipmap)

import sys
import clipshape_core_gene as csc

#functions

#maps clip peak ranges to genes
#needs to change to operate on gene ranges
def clipMapRanges(gene_range,clipdict):
    clipped_dict={}
    clipchr=gene_range[1]
    #gene_id,chr,strand,start,end
    for pos in range(gene_range[3],gene_range[4]):
        isclip=0
        try:
            for rng in clipdict[clipchr]:
                if pos >= rng[0] and pos <= rng[1]:
                    isclip=1
                    break
            if isclip==1:
                clipped_dict[pos]=1
            else:
                clipped_dict[pos]=0
        except KeyError:
            #no clip data for chr
            sys.stderr.write('WARNING: %s is mapped to %s, which has no CLIP data available. Skipping...\n' % (gene_range[0],gene_range[1]))
    return gene_range + [clipped_dict]

#grabs shape data for a stitched transcript.  takes output of clipMap
#def clipMerge(clipped_gene,shape_list):
#    for shape in shape_list:
#        if shape[0] == clipped_gene[0]:
#            shapevals=shape[3:]
#            if len(shapevals) != clipped_gene[4]-clipped_gene[3]+1:
#                sys.stderr.write('ERROR:  length mismatch in %s. Pos: %s Shape: %s.\n' % (shape[0],clipped_gene[4]-clipped_gene[3],len(shapevals)))
#                return 'FAIL'
#            else:
#                if clipped_gene[2] == '+':
#                    return clipped_gene + [shapevals]
#                elif clipped_gene[2] == '-':
#                    return clipped_gene + [shapevals[::-1]]
#            break

def clipMerge(clipped_gene,shape_data):
    try:
        shapevals=shape_data[clipped_gene[0]]
    except KeyError:
        sys.stderr.write('ERROR:  no shape data for %s. Pos: %s Shape: %s.\n' % (clipped_gene[0],clipped_gene[4]-clipped_gene[3],len(shapevals)))
        return 'FAIL'
    if len(shapevals) != clipped_gene[4]-clipped_gene[3]+1:
        sys.stderr.write('ERROR:  length mismatch in %s. Pos: %s Shape: %s.\n' % (shape[0],clipped_gene[4]-clipped_gene[3],len(shapevals)))
        return 'FAIL'
    else:
        if clipped_gene[2] == '+':
            return clipped_gene + [shapevals]
        elif clipped_gene[2] == '-':
            return clipped_gene + [shapevals[::-1]]

#final merge
def mergeAll(shaped_gene):
    finallist=[]
    clipped_pos=shaped_gene[5].keys()
    #snips off the noninclusive end
    positions=range(shaped_gene[3],shaped_gene[4])
    for i in range(len(positions)):
        tmplist=[shaped_gene[0], shaped_gene[1], shaped_gene[2], positions[i], shaped_gene[6][i]]
        if tmplist[3] in clipped_pos:
            tmplist.append(shaped_gene[5][tmplist[3]])
        else:
            tmplist.append('0')
        finallist.append(tmplist)
    return finallist

#catch chrs with no clip data
def toClip(chr_name,chr_dict):
    if chr_name in chr_dict.keys():
        return chr_dict[chr_name]
    else:
        return 'noclip'


###############################################
#inputs

shapefile=sys.argv[1]
clipfile=sys.argv[2]
gtffile=sys.argv[3]
outprefix=sys.argv[4]

#import gtf
sys.stderr.write('INFO: importing gtf...\n')
gtfdict=csc.importGTFGene(gtffile)

#import clip data
#this works for bedfiles
#or any 0-based tab-del format starting with chr,start,end
sys.stderr.write('INFO: importing clip data...\n')
clip_data={}
with open(clipfile) as infile:
    for line in infile:
        cliptmp=line.strip().split('\t')
        try:
            clip_data[cliptmp[0]].append((int(cliptmp[1])+1,int(cliptmp[2])+1))
        except KeyError:
            clip_data[cliptmp[0]]=[(int(cliptmp[1])+1,int(cliptmp[2])+1)]

#import shape data
sys.stderr.write('INFO: importing SHAPE data...\n')
#shape_data=[]
#with open(shapefile) as infile:
#    for line in infile:
#        shape_data.append(line.strip().split('\t'))

shape_data={}
with open(shapefile) as infile:
    for line in infile:
        shapetmp=line.strip().split('\t')
        shape_data[shapetmp[0]]=shapetmp[3:]

############################################


ranges_outname = "%s.generanges" % outprefix
f1=open(ranges_outname,'w')

clip_outname = "%s.clipmap_gene" % outprefix
f2=open(clip_outname,'w')

#id_list=[x[0] for x in shape_data]
id_list=shape_data.keys()

sys.stderr.write('INFO: Processing genes...\n')
genct=0
for id in id_list:
    #no longer necessary; these are genes
    #trans_stitched=csc.stitchTrans(id,gtfdict)
    #also not necessary
    #trans_ranges=csc.outputRanges(trans_stitched)
    #but I do need to reformat the gene range
    gene_range=csc.outputGeneRanges(id,gtfdict)

    gene_clipped=clipMapRanges(gene_range,clip_data)
    if gene_clipped != 'FAIL':
        gene_shaped=clipMerge(gene_clipped,shape_data)
        if gene_shaped != 'FAIL':
            gene_merged=mergeAll(gene_shaped)

            gene_out='\t'.join([str(x) for x in gene_range])
            f1.write(gene_out+'\n')

            clip_out=csc.flattenList(gene_merged)
            f2.write(clip_out+'\n')
            genct+=1
            if genct % 100 == 0:
                sys.stderr.write('INFO: %s genes processed successfully...\n' % str(genct))

f1.close()
f2.close()
