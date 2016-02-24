"""
sam2ng.py
Zhipeng Lu, 2015-10-16

Read in sam file with DG (duplex group) and XG (chiastic) tags,
assemble the DG into NG (non-overlapping groups).

testing a new version that minimize vertical space usage.
This version might be a little slower than the normal version using greedy algorithm.


test command:
python ~/Documents/scripts/paris/sam2ng.py AMT_Stress_trim_nodup_bc07_starhsXIST_l15p2.support_sorted.sam AMT_Stress_trim_nodup_bc07_starhsXIST_l15p2.support_sorted_NG.sam

"""

import sys, re
if len(sys.argv) < 3:
    print "sam2ng.py: assembles NG sam file from DG"
    print "also removes DGs with only duplicate reads and XG:i:2."
    print "Usage: python sam2ng.py inputsam outputsam"
    sys.exit()

inputsam = sys.argv[1]
outputsam = sys.argv[2]
inputsamh = open(inputsam, 'r')
outputsamh = open(outputsam, 'w')

samheader = ''
numreads = 0
dgdict = {} #dgname: [dginfo, dgreads ...], dginfo is a list of three elements
for line in inputsamh: #construct a dictionary with all DGs
    if line[0] == "@":
        samheader += line
        continue
    numreads += 1
    record = line.strip('\n').split()
    readstart = int(record[3])
    cigar = record[5]
    dgname = record[20]
    cigarbits = re.findall('\d+[MIN]', cigar)
    readspan = 0
    md = record[16].split(":")[-1]
    mdnum = len(re.findall('\d+[ATCG]', md))
    if mdnum > 1: continue #remove reads with >1 mismatches
    for i in cigarbits:
        if i[-1] == "I": readspan -= int(i[0:-1])
        else: readspan += int(i[0:-1])
    readend = readstart + readspan
    if dgname not in dgdict.keys():
        dginfo = [readstart, readend, 1]
        dgdict[dgname] = [dginfo, line]
    else:
        dgdict[dgname][0][0] = min(dgdict[dgname][0][0], readstart) #update the dgspan
        dgdict[dgname][0][1] = max(dgdict[dgname][0][1], readend)
        dgdict[dgname][0][2] += 1
        dgdict[dgname].append(line)

#for dgname in dgdict.keys(): print dgdict[dgname][0], len(dgdict[dgname])
#remove DG where all the reads have identical breaks in the CIGAR strings
#simply compare the N substrings for now
dgnamelist = dgdict.keys()
for dgname in dgnamelist:
    dgreads = dgdict[dgname][1:]
    breaklist = []
    for read in dgreads:
        record = read.strip('\n').split()
        cigar = record[5]
        cigarbits = tuple(re.findall('\d+[N]', cigar))
        breaklist.append(cigarbits)
    if len(list(set(breaklist))) == 1 : dgdict.pop(dgname)
    
outputsamh.write(samheader)
print "\nNumber of reads:", numreads
print "Number of duplex groups:", len(dgdict), "\nStarting NG assembly ...\n"



def addnginfo(dg, ngname): #add NG information (NG:i:x) to each read in the DG
    for k in range(1, len(dg)):
        #print len(dg)
        dg[k] = dg[k].strip('\n') + "\t" + ngname
    return dg


dglist = dgdict.values() #sort dgdict to a list based on number of reads in each DG.
dglistsorted = [dg for dg in sorted(dglist, key=lambda dg: dg[0][2], reverse=True)]
ng0 = dglistsorted.pop(0)
ng0 = [[ng0[0]]] + ng0[1:]
ngname = "NG:i:1" #start from 1
ng0 = addnginfo(ng0, ngname)
nglist = [ng0] #example: [[[1,20],read1,read2], [[10,40], read3,read4]]
ngindex = 2 #track and name NG
#print ng0[0]
#print len(ng0)

"""
#This part uses greedy algorithm to pack NG, fast but not efficient in vertical space
for dg in dglistsorted:
    for ng in nglist:           #dg/ng: [[1,20],read1,read2]
        if dg[0][0] > ng[0][1] or dg[0][1] < ng[0][0]:
            ng[0][0] = min(dg[0][0], ng[0][0])
            ng[0][1] = max(dg[0][1], ng[0][1])
            ngname = ng[1].split('\t')[-1]
            dg = addnginfo(dg, ngname)
            ng += dg[1:]
            break
    else:
        ngname = "NG:i:" + str(ngindex)
        dg = addnginfo(dg, ngname)
        nglist.append(dg) #start a new NG with the DG
        ngindex +=1
print "Number of assembled non-overlapping groups:", ngindex
#####################################################################################
"""

# ng: [[[1,20], [50, 60]],read1,read2]
# dg: [[10,30],read5,read6]
def overlap(dg, ng):
    for interval in ng[0]: 
        if interval[0] <= dg[0][0] <= interval[1] or \
           interval[0] <= dg[0][1] <= interval[1] or \
           dg[0][0] <= interval[0] <= dg[0][1] or \
           dg[0][0] <= interval[1] <= dg[0][1] :
            return 1
            break
    else: return 0


#assemble the DG to NG
for dg in dglistsorted:
    for ng in nglist:           #dg/ng: [[1,20],read1,read2]
        if not overlap(dg, ng): 
            ng[0].append(dg[0])
            ngname = ng[1].split('\t')[-1]
            dg = addnginfo(dg, ngname)
            ng += dg[1:]
            break
    else:
        ngname = "NG:i:" + str(ngindex)
        dg = addnginfo(dg, ngname)
        newng = [[dg[0]]] + dg[1:]
        nglist.append(newng) #start a new NG with the DG
        ngindex +=1
        
print "Number of assembled non-overlapping groups:", ngindex




allng = ''
for ng in nglist:
    #print "Number of reads in", ng[1].split()[-1], "is:\t", len(ng)
    ngout = '\n'.join(ng[1:])
    allng += (ngout + '\n')
outputsamh.write(allng)
      

inputsamh.close()
outputsamh.close()


        
"""
samtools view -bS -o a.bam MALAT1_AMT_Stress_trim_nodup_norm_starhg38_l15p2_geometric_common25_NGmin.sam
samtools sort a.bam a_sorted
samtools index a_sorted.bam

"""
