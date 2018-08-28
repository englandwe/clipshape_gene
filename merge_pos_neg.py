
import sys

posfile=sys.argv[1]
negfile=sys.argv[2]

with open(posfile) as infile:                                                                   #no     #yes    #no
#ENST00000307046        12      -       102417817       NULL    1       86      A       81      -43     -54     AGCAA
   for line in infile:
#       print(line)
       tmpline=line.strip().split('\t')
       print('\t'.join(tmpline[0:9] + [tmpline[10],'pos']))


with open(negfile) as infile:                                                                 #new col: [6] -7
#ENST00000491035        X       +       154400128       NULL    0       0       G       0
#ENST00000585635	18	-	36828683	36829036	180	195	1	T	0.013
#0,1,2,3+7-int(1),9,'0',7,8,?,7-8
    for line in infile:
        tmpline=line.strip().split('\t')
        outline=tmpline[0:3]+[int(tmpline[3])+int(tmpline[7])-1,tmpline[9],'0',tmpline[7],tmpline[8],'?',int(tmpline[7])-253,'neg']
        print('\t'.join([str(x) for x in outline]))

