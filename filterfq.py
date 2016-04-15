# this program filters fastq files to remove short reads

import sys


infile  = sys.argv[1]
outfile = sys.argv[2]

infileh = open(infile, 'r')
outfileh= open(outfile, 'w')


while True:
    lines = []
    for i in range(4):
        try:
            lines.append(infileh.next())
        except StopIteration:
            infileh.close()
    if len(lines[1].strip('\n')) > 30:
        outfileh.write(lines)
outfileh.close()
