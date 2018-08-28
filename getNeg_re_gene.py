#takes: 
#list of tab-separated positive ranges (chr, start, stop, strand). start should always be smaller than stop.  stop is not inclusive.
#genomic fasta
#txranges file (from clipmap).  genomic ranges for all exons. format: tx chr strand start stop
#shape-like output file

#returns: shape-like data for an equal number of negative control ranges


import sys
#import os
import random
import re
from Bio import SeqIO,motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
#from itertools import groupby
#from operator import itemgetter
import clipshape_core_gene as csc

def makeTxRanges(gtffile,shapedict):
    tx_ranges=[]
    sys.stderr.write('Importing GTF...\n')
    gtfdict=csc.importGTF(gtffile)
    sys.stderr.write('GTF imported. building tx set...\n')
    #tx_set=set([x.attdict['transcript_id'] for x in gtf_list if x.feature == 'transcript'])
    #sys.stderr.write('tx set built.  length:%s...\n' % len(tx_set))
    txct=1
    for tx in shapedict.keys():
        sys.stderr.write('Processing tx %s...\n' % txct)
        trans_stitched=csc.stitchTrans(tx,gtfdict)
        sys.stderr.write('Stitched tx %s...\n' % txct)
        trans_ranges=csc.outputRanges(trans_stitched)
        sys.stderr.write('Got ranges for tx %s...\n' % txct)
        tx_ranges+=trans_ranges
        txct+=1
    sys.stderr.write('ranges made!\n')    
    return tx_ranges

#import positive ranges
def importPosRanges(posfile):
    rnglist=[]
    with open(posfile) as infile:
        for line in infile:
            #chr, start, stop (not inclusive), strand
            posrng=line.strip().split('\t')
            rnglist.append([posrng[0],int(posrng[1]),int(posrng[2]),posrng[3]])
    return rnglist

#import known tx ranges.  useful if you've already run clipmap.
def importTxRanges(txfile):
    tx_ranges=[]
    with open(txfile) as inrange:
        for line in inrange:
            tx_ranges.append(line.strip().split('\t'))
    #tx,chr,strand,start,stop
    return tx_ranges

def importShape(shapefile):
    shapedict={}
    with open(shapefile) as infile:
        for line in infile:
            tmpline=line.strip().split('\t')
            shapedict[tmpline[0]]=tmpline[3:]
    return shapedict

#get sequence data for a rangelist
def getSeq(rnglist,fasta_data):
    new_rnglist=list(rnglist)
    for rng in new_rnglist:
        #print(rng)
        if rng[3] == '+':
            #find that chromosome(it'll be the dict ID in fasta dict) and get the sequence for those positions
            range_seq=fasta_data[rng[0]].seq[rng[1]-1:rng[2]-1]
        elif rng[3] == '-':
            range_seq=fasta_data[rng[0]].seq[rng[1]-1:rng[2]-1].reverse_complement()
        rng.append(range_seq)
        #print(range_seq)
    #chr, start, stop (not inclusive), strand, seq    
    return new_rnglist

#make a motif from the positive ranges    
def getMotif(new_rnglist):
    motif_list=[]
    for entry in new_rnglist:
        motif_seq=entry[4]
        motif_list.append(motif_seq)
    #make motif.  make lc uppercase
    clip_motif=motifs.create([x.upper() for x in motif_list if 'N' not in x],alphabet=IUPAC.unambiguous_dna)
    #set pseudocounts
    clip_motif.pseudocounts=0.5
    return clip_motif
    
def findUnbound(tx_ranges,clip_motif,fasta_data,clip_ranges,shapedict,fp_offset,tp_offset,peakcount):
    true_neg=0
    negout=[]
    if peakcount=='auto':
        peakcount=len(clip_ranges)
        sys.stderr.write('%\n' % peakcount)
    else:
	peakcount=int(peakcount)
    #clip_pssm=clip_motif.pssm
    motif_len=len(re.sub('\[[ACTG]+\]','N',clip_motif))
    random.shuffle(tx_ranges)
    txct=0
    #tx,chr,strand,start,stop (inclusive)
    #for each tx range (aka exon), get sequence
    for txrng in tx_ranges:
        txct+=1
        if txct%100 == 0:
            sys.stderr.write('testing range %s\n' % (txct))
        if txrng[2] == '+':
            range_seq=fasta_data[txrng[1]].seq[int(txrng[3])-1:int(txrng[4])]
        elif txrng[2] == '-':
            range_seq=fasta_data[txrng[1]].seq[int(txrng[3])-1:int(txrng[4])].reverse_complement()
        #see if there are any hits
        #sys.stderr.write(str(range_seq)+'\n')
        tmpsites=[]
        #for position,score in clip_pssm.search(range_seq,threshold=3.0):
            #tmpsites.append([position,score])
        for m in re.finditer(clip_motif,str(range_seq)):
            tmpsites.append([m.start(),m.group()])
        #pointless to go further if there are no appropriate hits
        tmpsites_pos=[x for x in tmpsites if x[0] >= 1]
        #if there are hits, great!  grab an offset around the hit and add it to the list
        #hold on though; we only want one per exon, to spread out the hits.  so pick one at random
        if len(tmpsites_pos) > 0:
            random.shuffle(tmpsites_pos)
            chosen_one=tmpsites_pos[0]
            negstart=chosen_one[0]-fp_offset
            negend=chosen_one[0]+tp_offset+motif_len
            #make sure you can get the entire requested range
            if negstart >= 1 and negend <= len(range_seq):
                #final check:  does this range overlap a positive range
                #get genomic coords
                if txrng[2] == '+':
                    motifstart_genomic=int(txrng[3])+negstart-1
                    motifstop_genomic=int(txrng[3])+negend
                elif txrng[2] == '-':
                    #get distance of end of motif from right end of string
                    motifstop_genomic=int(txrng[4])-negstart
                    motifstart_genomic=int(txrng[4])-negend-1
                    #motifstart_genomic=len(range_seq)-negstart-len(clip_motif)-offset
                    #motifstop_genomic=motifstart_genomic+len(clip_motif)+offset+offset
                sys.stderr.write('motif range:%s,%s\n' % (motifstart_genomic,motifstop_genomic))
                #check for range overlaps
                isPeak=0
                for clip_peak in clip_ranges:
                    #chr, start, stop (not inclusive), strand
                    for pos in range(motifstart_genomic,motifstop_genomic):
                        if pos >=clip_peak[1] and pos < clip_peak[2]:
                            isPeak=1
                if isPeak == 1:
                    sys.stderr.write('that was a peak\n')
                else:
                    ##
                #if clipped_gene[2] == '+':
                #    return clipped_gene + [shapevals]
                #elif clipped_gene[2] == '-':
                #    return clipped_gene + [shapevals[::-1]]
                    sys.stderr.write('not a peak,excellent\n')
                    negseq=range_seq[negstart:negend].upper()
                    if txrng[2] == '+':
                        shapevals=shapedict[txrng[0]][negstart:negend]
                    elif txrng[2] == '-':
                        shapevals=shapedict[txrng[0]][::-1][negstart:negend]
                    negout.append(txrng+[negstart,negend,negseq,shapevals])
                    true_neg+=1
                    sys.stderr.write('found a negative range, running total is %s\n' % (true_neg))
        if true_neg == peakcount:
            break
    return negout

#['ENST00000612958', 'X', '+', '49338836', '49338952', 86, 109, \
#Seq('ccaataaagctttacagccttct', IUPACUnambiguousDNA()), \
#['0.188', '0.014', '0.127', '0.351', '0.029', '0.062', '0.223', ....]
#tx,chr,strand,start,stop,txrngstart,txrngstop,seq,shape
def outfmt(outlist):
    final_out=[]
    for item in outlist:
        posct=1
        for idx,shape in enumerate(item[8]):
            #outitem=item[0:7]+[posct,str(item[7])[idx],item[8][idx]]
            outitem=item[0:3]+[int(item[3])+idx,'blank',int(item[5]+idx),'blank']+[posct,str(item[7])[idx],item[8][idx]]
            final_out.append(outitem)
            posct+=1
    return final_out

#posfile='polyA_ranges_35up_shape'
#txfile='auto'
#gtffile='../../../ref_txomes/human_rRNA/Homo_sapiens.GRCh38.88.plusrRNA.gtf'
#fastafile='../../../ref_txomes/human_rRNA/Homo_sapiens.GRCh38.dna_sm.primary_assembly.plusrRNA.fa'
#shapefile='icshape.unfilt.out'
#peakcount='auto'
#fp_offset=101-36
#tp_offset=0

#filename or the actual motif
posfile=sys.argv[1]
#auto or motif
motifseq=sys.argv[2]
#auto or filename
txfile=sys.argv[3]
#only needed if txfile=='auto'
gtffile=sys.argv[4]
fastafile=sys.argv[5]
shapefile=sys.argv[6]
#integer or auto
peakcount=sys.argv[7]
fp_offset=int(sys.argv[8])
tp_offset=int(sys.argv[9])

motif_alphabet=set('ACGT[]')

sys.stderr.write('importing fasta\n')
fasta_dict=SeqIO.to_dict(SeqIO.parse(fastafile, "fasta", alphabet=IUPAC.unambiguous_dna))
sys.stderr.write('importing SHAPE\n')
shapedict=importShape(shapefile)

if txfile == 'auto':
    sys.stderr.write('generating txranges\n')
    #generate transcript ranges
    tx_ranges=makeTxRanges(gtffile,shapedict)
else:
    sys.stderr.write('importing txranges\n')
    tx_ranges=importTxRanges(txfile)

sys.stderr.write('importing positive ranges\n')
pos_ranges=importPosRanges(posfile)

if motifseq == 'auto':
    #Get positive seq
    sys.stderr.write('getting positive seqs\n')
    pos_ranges_withseq=getSeq(pos_ranges,fasta_dict)
    #Motif
    sys.stderr.write('finding motif\n')
    pos_motif=getMotif(pos_ranges_withseq)
    #sys.stderr.write(str(len(pos_motif))+'\n')
    #sys.stderr.write(str(pos_motif.pwm)+'\n')
    #sys.stderr.write(str(pos_motif.pssm)+'\n')
elif set(motifseq).issubset(motif_alphabet):
    #use that as the motif directly
    pos_motif=motifseq
else:
    sys.stderr.write('thats not auto or ACTG[]\n') 

#motif_list=['TTT','TTT','TTT','TTT','TTT','TTT','TTT','TTT','TTT','TTT']
#pos_motif=motifs.create([x.upper() for x in motif_list if 'N' not in x],alphabet=IUPAC.unambiguous_dna)
#pos_motif.pseudocounts=0.5
unbound=findUnbound(tx_ranges,pos_motif,fasta_dict,pos_ranges,shapedict,fp_offset,tp_offset,peakcount)
print(csc.flattenList(outfmt(unbound)))
#print(outlist)




