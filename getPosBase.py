import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

clipmap_file=sys.argv[1]
fastafile=sys.argv[2]

sys.stderr.write('importing fasta...\n')
fasta_dict=SeqIO.to_dict(SeqIO.parse(fastafile, "fasta", alphabet=IUPAC.unambiguous_dna))

rcdict={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','g':'c','c':'g'}
def getBase(fasta_data,rcdict,chr,strand,pos):
    if strand == '+':
        base=fasta_data[chr].seq[pos-1]
    elif strand == '-':
        base=rcdict[fasta_data[chr].seq[pos-1]]
    return base

#ENST00000425361        2       +       130190577       0.704   NEWVAL	1
sys.stderr.write('commencing base addition...\n')
with open(clipmap_file) as infile:
    peakcounter=0
    poscounter=0
    prevline=next(infile).strip().split('\t')
    if prevline[-1] == '1':    
        posbase=getBase(fasta_dict,rcdict,prevline[1],prevline[2],int(prevline[3]))
        outline=prevline+[poscounter,posbase,peakcounter]
        print('\t'.join([str(x) for x in outline]))
        poscounter+=1
    for line in infile:
        #sys.stderr.write('peakcount: %s\n' % peakcounter)
        tmpline=line.strip().split('\t')
        if tmpline[-1] == '1':
            if tmpline[0:3] == prevline[0:3] and int(tmpline[3]) == int(prevline[3])+1:
                #consecutive
                posbase=getBase(fasta_dict,rcdict,tmpline[1],tmpline[2],int(tmpline[3]))
                outline=tmpline+[poscounter,posbase,peakcounter]
                print('\t'.join([str(x) for x in outline]))
                prevline=list(tmpline)
                poscounter+=1
            else:
                poscounter=0
                peakcounter +=1
                posbase=getBase(fasta_dict,rcdict,tmpline[1],tmpline[2],int(tmpline[3]))
                outline=tmpline+[poscounter,posbase,peakcounter]
                print('\t'.join([str(x) for x in outline]))
                prevline=list(tmpline)
                poscounter+=1
    
