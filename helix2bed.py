"""
This script converts manual helix annotations to bed formats

python helix2bed.py file.helix file.bed chrom adjust rnasize rnastrand

The adjust parameter describes the offset between the input and output (input-output)

Helix annotations (6 fields, 1-based): 
Helix	1start	1end	2start	2end	UU

Bed format (12 fields, 0-based): 
chrom	chromStart	chromEnd	name	score   strand	thickStart
thickEnd	itemRgb	blockCount	blockSizes	blockStarts
"""
import sys
print "Usage: python helix2bed.py file.helix file.bed chrom adjust rnasize rnastrand"

helixfile = sys.argv[1]
bedfile = sys.argv[2]
chrom = sys.argv[3]
adjust = int(sys.argv[4])
rnasize = int(sys.argv[5])
rnastrand = sys.argv[6]

helixf = open(helixfile, 'r')
bedf = open(bedfile, 'w')
beds = []


for line in helixf:
    helix = line.strip('\n').split()
    name = helix[0]
    strand = rnastrand
    p1, p2, p3, p4 = int(helix[1]), int(helix[2]), int(helix[3]), int(helix[4])
    if strand == '-':
        p1, p2, p3, p4 = rnasize - p4, rnasize - p3, rnasize - p2, rnasize - p1
        
    chromStart = str(p1 + adjust)
    chromEnd = str(p4 + adjust)
    score = '1000'
    thickStart = chromStart
    thickEnd = chromStart
    itemRgb = '0,0,0'
    blockCount = '2'
    blockSizes = str(p2 - p1 + 1) + ',' + str(p4 - p3 + 1)
    blockStarts = '0,' + str(p3 - p1)
    bed = [chrom, chromStart, chromEnd, name, score, "*", \
           thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts]
    beds.append(bed)

 
for bed in beds:
    bedf.write("\t".join(bed)+"\n")
    

helixf.close()
bedf.close()

