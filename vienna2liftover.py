"""""""""
This program converts Clustal alignments in vienna format
to liftOver chain format.

Zhipeng Lu, 2015-08-23

python vienna2liftover.py aln.vienna aln.liftoverchain

"""""""""

import sys
if len(sys.argv) < 3: 
    print "Usage: python vienna2liftover.py aln.vienna aln.liftoverchain"
    sys.exit()

alnvienna = sys.argv[1]
alnliftover = sys.argv[2]
alnviennaf = open(alnvienna, 'r')
alnliftoverf = open(alnliftover, 'w')

aln1name = alnviennaf.next()[1:-1]
aln1seq  = alnviennaf.next().strip('\n')
aln2name = alnviennaf.next()[1:-1]
aln2seq  = alnviennaf.next().strip('\n')
alnviennaf.close()
seqlen = len(aln1seq)
aln1size = seqlen - aln1seq.count('-')
aln2size = seqlen - aln2seq.count('-')

outlist = ''
block = [0,0,0]
#convert nucleotides to N
nadict = {'A':'N', 'C':'N', 'G':'N', 'T':'N', 'U':'N', 'N':'N', '-':'-'}


#alignment blocks always start from 'NN' and end with 'NN'
i = tlead = qlead = 0
currentpos = nadict[aln1seq[0]] + nadict[aln2seq[0]]
while currentpos != 'NN':
    if currentpos == 'N-': tlead += 1
    else: qlead += 1
    i += 1
    currentpos = nadict[aln1seq[i]] + nadict[aln2seq[i]]


#process from 'NN' till the end
prevpos = 'NN'
while i < seqlen:
    currentpos = nadict[aln1seq[i]] + nadict[aln2seq[i]]
    i += 1
    if (prevpos != 'NN') and (currentpos == 'NN'):
        outlist += (str(block[0]) + "\t" + str(block[1]) + "\t" + str(block[2]) + "\n")
        block = [1,0,0]
    else:
        if currentpos == 'NN':   block[0] += 1
        elif currentpos == 'N-': block[1] += 1
        else: block[2] += 1
    prevpos = currentpos   
lastmatch = block[0]
ttrail    = block[1]
qtrail    = block[2]

    
# make the following header line
# chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
headers = ['chain', '10000', aln1name, str(aln1size), '+', str(tlead), str(aln1size-ttrail), \
           aln2name, str(aln2size), '+', str(qlead), str(aln2size-qtrail), '1']
header = "\t".join(headers)
alnliftoverf.write(header + "\n")
alnliftoverf.write(outlist + str(lastmatch) + '\n')
alnliftoverf.close()


"""
>seq1
NNNN--NNNN--N--
>seq2
--NNNNNN--NNNNN

len = 15
len(seq1) = 9
len(seq2) = 11
"""
