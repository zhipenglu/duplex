"""
liftover_clean.py

This script removes the duplicate coordinates from liftover output.
The liftover output should be sorted by coordinates first. 

python ~/Documents/scripts/liftover_clean.py icSHAPE_PAinvivo_hg38_sorted.bedgraph \
icSHAPE_PAinvivo_hg38_sorted_clean.bedgraph
bedGraphToBigWig icSHAPE_PAinvivo_hg38_sorted_clean.bedgraph \
~/Documents/chang/hg38/hg38_size.txt icSHAPE_PAinvivo_hg38_sorted_clean.bigwig
"""

import sys

if len(sys.argv) < 3:
    print "Usage: python liftover_clean.py bedgraphin bedgraphout"
    sys.exit()

bedgraphin  = open(sys.argv[1], 'r')
bedgraphout = open(sys.argv[2], 'w')

lastend = 0
lastchrom = ''
outstring = ''

for line in bedgraphin:
    interval = line.split()[0:3]
    chrom = interval[0]
    start = int(interval[1])
    if chrom != lastchrom or start >= lastend:        
        lastend = int(interval[2])
        lastchrom = chrom
        outstring += line

bedgraphout.write(outstring)
bedgraphin.close()
bedgraphout.close()
