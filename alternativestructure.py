"""
this script go through the entire bed file of DGs
and extract potential pseudoknots (interlocking DGs)
and alternative structures

requires RNAcofold in proper path:
for my own Mac: RNAcofold
for Changrila: /home/zhipeng/bin/ViennaRNA/bin/RNAcofold

install intervaltree on changrila:
pip install --install-option="--prefix=/home/zhipeng/lib" intervaltree

test command:
cd /Users/lu/Documents/chang/psoralen/higherorder/
python ~/Documents/scripts/duplex/alternativestructure1.py \
AMT_Stress_trim_nodup_starhsXIST_l15p2_NGmin1386_fix.bed \
/Users/lu/Documents/chang/psoralen/examples/XIST/hsXIST_igv/hsXIST.fa \
AMT_Stress_trim_nodup_starhsXIST_l15p2_NGmin1386_fix.alt

Number of all DGs: 1386
Length of all intervals: 15599
Number of armpairs: 13609
Number of potential alternative pairs: 3462
cut -f4,17 AMT_Stress_trim_nodup_starhsXIST_l15p2_NGmin1386_fix.alt \
| tr '\t' '\n' | sort -u | wc -l
Out of the 1386 DGs, 1252 of them are involved in alternative structures
some of them are involved in up to 20 alternative structure pairs. 



a shorter example
cd /Users/lu/Documents/chang/psoralen/higherorder/
python ~/Documents/scripts/duplex/alternativestructure.py a.bed \
/Users/lu/Documents/chang/psoralen/examples/XIST/hsXIST_igv/hsXIST.fa alt

For transcriptomic analysis on Changrila (analyze all 8 samples separately):
cd /home/zhipeng/HEK1/HEK1_hg38/HEK1_bc07/
python ~/bin/alternativestructure.py \
AMT_Stress_trim_nodup_bc07_starhg38_geometricfilterfix.bed /home/zhipeng/hg38/hg38.fa \
AMT_Stress_trim_nodup_bc07_starhg38_geometricfilterfix.alt > \
AMT_Stress_trim_nodup_bc07_starhg38_geometricfilterfix.altlog &

cd /home/zhipeng/HeLa6/HeLa6_hg38/
python ~/bin/alternativestructure.py \
AMT_HeLa6_trim_nodup_starhg38_geometricfilterfix.bed /home/zhipeng/hg38/hg38.fa \
AMT_HeLa6_trim_nodup_starhg38_geometricfilterfix.alt > \
AMT_HeLa6_trim_nodup_starhg38_geometricfilterfix.altlog &

cd /home/zhipeng/mES4/Mettl3_mm10/
python ~/bin/alternativestructure.py \
AMT_Mettl3_trim_nodup_bc01_starmm10_geometricfilterfix.bed \
/home/zhipeng/mm10/mm10upper.fa \
AMT_Mettl3_trim_nodup_bc01_starmm10_geometricfilterfix.alt > \
AMT_Mettl3_trim_nodup_bc01_starmm10_geometricfilterfix.altlog &


"""

import sys, itertools, subprocess, datetime
sys.path.append("/home/zhipeng/lib/lib/python2.7/site-packages") # for use on changrila
import intervaltree

if len(sys.argv) < 4:
    print "Usage: python alternativestructure.py bedfile fastafile alt_out"
    print "Overlapping thresholds are set in this script."
    sys.exit()

#TO DO: remove the size limit on alternative structures (pksize). 
#TO DO: read in genome fasta files for look up of the sequence.
#TO DO: test conflicts of base pairing with simple bp predictions

pksize = 1000000000   
over_threshold = 0.5 #overlap of overlapped arm in alternative structures
alt_threshold = 0.2 #overlap of the other arm in alternative structures
coveragemin = 2 #threshold number of reads for each DG
bedfile = open(sys.argv[1], 'r')
fastafile = open(sys.argv[2], 'r')
altoutputfile = open(sys.argv[3], 'w')

########read the fasta files to test 
def readfasta(fastafile):
    fastadict = {}
    fh = fastafile
    faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        fastadict[header] = seq
    return fastadict
fastadict = readfasta(fastafile)
print "\n", str(datetime.datetime.today())
print "Number of references:", len(fastadict)

########construct a DG tree for each chrom, the bed file may not be sorted
dgtrees = {} #each dgtree is named after the chrom name
intervals = {} #positions for each chrom is stored in a list
for line in bedfile:
    bed = line.strip('\n').split()
    blockSizes = bed[10].strip(',').split(',')
    chrom = bed[0]
    begin1, end1 = int(bed[1]), int(bed[1]) + int(blockSizes[0])
    begin2, end2 = int(bed[2]) - int(blockSizes[1]), int(bed[2]) 
    if chrom not in dgtrees:
        dgtrees[chrom] = intervaltree.IntervalTree()
        intervals[chrom] = []
    dgtrees[chrom].addi(begin1, end1, line)
    dgtrees[chrom].addi(begin2, end2, line)
    for i in range(begin1, end1):
        intervals[chrom].append(i)
    for i in range(begin2, end2):
        intervals[chrom].append(i)


#########convert dgtrees to armpairs
treesize, intervallength, armpairscount = 0, 0,0
chroms = dgtrees.keys()
armpairs = {}
#example armpair: (Interval(454, 475, 'bedrecord'), Interval(454, 476, 'bedrecord'))
for chrom in chroms:
    treesize += len(dgtrees[chrom])
    intervals[chrom] = sorted(list(set(intervals[chrom])))
    intervallength += len(intervals[chrom])
    #make all dgclusters from this chromosome
    dgclusters = set()
    for i in intervals[chrom]:
        dgcluster = tuple(dgtrees[chrom][i])
        if dgcluster and len(dgcluster) > 1: dgclusters.add(dgcluster)
    #make all armpairs from this chromosome
    armpairs[chrom] = set()
    for dgcluster in dgclusters:
        for dg1, dg2 in list(itertools.combinations(dgcluster, 2)):
            armpairs[chrom].add(tuple(sorted((dg1, dg2))))
    armpairscount += len(armpairs[chrom])


    
print "\n", str(datetime.datetime.today())
print "Number of all DGs:", treesize/2
print "Length of all intervals:", intervallength
print "Number of armpairs:", armpairscount

def overlap(al, ar, bl, br):
    #two intervals a [al, ar] and b [bl, br]
    #returns the ratio of overlap
    if ar < bl or al > br: return 0
    else:
        alength, blength = ar - al, br - bl
        coords = sorted([al, ar, bl, br])
        overlaplength = float(coords[2] - coords[1])
        overlap = min(overlaplength/alength, overlaplength/blength)
        return overlap
        
def alternativecheck(intvl1, intvl2, over_threshold, alt_threshold):
    #returns whether the two intervals overlap
    #use the overlap() function defined above
    if (overlap(intvl1[0], intvl1[1], intvl2[0], intvl2[1]) > over_threshold) and \
       (overlap(intvl1[2], intvl1[3], intvl2[2], intvl2[3]) < alt_threshold) or \
       (overlap(intvl1[0], intvl1[1], intvl2[0], intvl2[1]) < alt_threshold) and \
       (overlap(intvl1[2], intvl1[3], intvl2[2], intvl2[3]) > over_threshold) or \
       (overlap(intvl1[0], intvl1[1], intvl2[2], intvl2[3]) > over_threshold) or \
       (overlap(intvl1[2], intvl1[3], intvl2[0], intvl2[1]) > over_threshold):
        return 1
    else: return 0

def bpcheck(chrom, intvl1, intvl2, fastadict):
    #use RNAcofold to obtain structure and test overlap
    #use the overlap() function defined above
    #return predicted structure
    structures = []
    si = [] #structure intervals: [[a1, a2, a3, a4], [b1, b2, b3, b4]]
    for intvl in [intvl1, intvl2]: 
        seq = fastadict[chrom][intvl[0]:intvl[1]]+"&"+fastadict[chrom][intvl[2]:intvl[3]]
        constraint = "<" * (intvl[1] - intvl[0]) + ">" * (intvl[3] - intvl[2])
        p = subprocess.Popen(["RNAcofold", "-C"], \
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        helix = p.communicate(input = seq + "\n" + constraint)[0].split('\n')
        structure = helix[1][0:len(helix[0])]
        structures.append(structure)
        structure = helix[1][0:len(helix[0])].replace(")", "(").split("&")
        structureintvl = [0,0,0,0]
        #remove the dots at the begining and end to make the new intervals for brackets. 
        if structure[0][0] != ".": structureintvl[0] = intvl[0]
        else: structureintvl[0] = intvl[0] + len(structure[0].split("(")[0])
        if structure[0][-1] != ".": structureintvl[1] = intvl[1]
        else: structureintvl[1] = intvl[1] - len(structure[0].split("(")[-1])
        if structure[1][0] != ".": structureintvl[2] = intvl[2]
        else: structureintvl[2] = intvl[2] + len(structure[1].split("(")[0])
        if structure[1][-1] != ".": structureintvl[3] = intvl[3]
        else: structureintvl[3] = intvl[3] - len(structure[1].split("(")[-1])
        si.append(structureintvl)
    maxoverlap = max(overlap(si[0][0], si[0][1], si[1][0], si[1][1]), \
           overlap(si[0][0], si[0][1], si[1][2], si[1][3]), \
           overlap(si[0][2], si[0][3], si[1][0], si[1][1]), \
           overlap(si[0][2], si[0][3], si[1][2], si[1][3]))
    #print maxoverlap, structures[0], structures[1]       
    if maxoverlap > over_threshold: return maxoverlap, structures[0], structures[1]
    else: return 0, 0, 0
    




########check each cluster of DG for potential pseudoknot and alternative structures
#hsXIST 840 1935 Group_56 3 * 840 840 0,0,0 2 20,20 0,1075
altpairs = []
altpairsbed = [] #bed12_1 + structure_1 + bed12_2 + structure_2 + maxoverlap.
for chrom in chroms:
    for armpair in armpairs[chrom]:
        dg1, dg2 = armpair # a tuple of two intervals
        dg1info = dg1.data.strip('\n').split("\t")
        dg2info = dg2.data.strip('\n').split("\t")
        start1 = int(dg1info[1])
        start2 = int(dg2info[1])
        end1 = int(dg1info[2])
        end2 = int(dg2info[2])
        blockSizes1 = [int(i) for i in dg1info[10].strip(',').split(',')]
        blockSizes2 = [int(i) for i in dg2info[10].strip(',').split(',')]
        blockStarts1 =[int(i) for i in dg1info[11].strip(',').split(',')]
        blockStarts2 =[int(i) for i in dg2info[11].strip(',').split(',')]
        intvl1 = [start1,start1+blockSizes1[0],end1-blockSizes1[1],end1]
        intvl2 = [start2,start2+blockSizes2[0],end2-blockSizes2[1],end2]
        ###test alternative structures
        maxoverlap, structure1, structure2 = bpcheck(chrom, intvl1, intvl2, fastadict)
        if alternativecheck(intvl1, intvl2, over_threshold, alt_threshold) and \
           min(int(dg1info[4]), int(dg2info[4])) >= coveragemin and maxoverlap: 
                altpairs.append(dg1info[3] + "\t" + dg2info[3])
                bedpair = [dg1.data.strip('\n'), structure1, \
                           dg2.data.strip('\n'), structure2, str(maxoverlap)]
                altpairsbed.append("\t".join(bedpair))

#for i in altpairs: print i
#for i in altpairsbed: print i
print "\n", str(datetime.datetime.today())
print "Number of potential alternative pairs:", len(altpairs)

altpairsout = "\n".join(altpairsbed) + "\n"
altoutputfile.write(altpairsout)
bedfile.close()
altoutputfile.close()

