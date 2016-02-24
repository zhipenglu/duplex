"""""""""""""""""""""
This program converts the ct file to bed file for display of RNA structures.
Zhipeng Lu
2015-08-20

Basic usage:
input:
input ct file  HS_LSU_3D_CT.ct
chr name       hs45S
input strand   "+"
RNA start      1
CT file format:
left	nucleotide	right-1	right	right+1	name

output:
bed file

current test command:
python ~/Documents/scripts/paris/ct2bed.py HS_LSU_3D_CT.ct hs45S "+" 1 rnasize test.bed


"""""""""""""""""""""

import sys
if len(sys.argv) < 7:
    print "Usage: python ct2bed.py file.ct chrname strand rnastart rnasize file.bed"
    sys.exit()
    
inputct   = sys.argv[1]
chrname   = sys.argv[2]
rnastrand = sys.argv[3]
rnastart  = int(sys.argv[4])
rnasize   = int(sys.argv[5])
outfile   = sys.argv[6]

outbeds = '' # make the output a single string for fast printing
inputctf = open(inputct, 'r')
outf = open(outfile, 'w')
outf.write("track graphType=arc\n")
headerline = next(inputctf) #skip first line, becareful about the first line
    
for line in inputctf:
    bp = line.split()
    if len(bp) < 5: continue
    bp1 = int(bp[0]) #one side of a basepair
    #print bp1
    bp2 = int(bp[4]) #the other side of the basepair
    if bp1 < bp2:
        if rnastrand == '-':
            bp1, bp2 = rnasize - bp2, rnasize - bp1
        outbeds += (chrname + "\t" + str(bp1+rnastart) + "\t" + \
                   str(bp2+rnastart) + "\t*\t1\t" + rnastrand + \
                    "\t" + str(bp1+rnastart) + "\t" + str(bp1+rnastart) + "\t0,0,0\n")
        
outf.write(outbeds)
inputctf.close()




