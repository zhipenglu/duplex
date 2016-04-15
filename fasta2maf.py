"""
This script converts the concatenated fasta back to maf in one block:
The coordinates start from 1 for every homolog.

"""

import sys

if len(sys.argv) < 3:
    print "Usage: python fasta2maf.py aligned.fasta aligned_cat.maf"
    sys.exit()

fastafile = sys.argv[1]
fastafileh = open(fastafile, 'r')
concatmaf = sys.argv[2]
concatmafh = open(concatmaf, 'w')


outrecord = '##maf version=1\na score=10000\n'
seqname = fastafileh.readline().strip('\n').strip('>')
seq = ''

for line in fastafileh:
    if line[0] == '>':
        outrecord += ("s " + seqname + '\t1\t' + str(len(seqclean)) + \
                      '\t+\t' + str(len(seqclean)) + '\t' + seq + '\n')
        seqname = line.strip('\n').strip('>')
        seq = ''
    else:
        seq += line.strip('\n')
        seqclean = seq.translate(None, '-')
outrecord += ("s " + seqname + '\t1\t' + str(len(seqclean)) + \
              '\t+\t' + str(len(seqclean)) + '\t' + seq + '\n')

concatmafh.write(outrecord)
fastafileh.close()
concatmafh.close()
