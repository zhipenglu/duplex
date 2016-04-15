"""
samfilterdg.py
Zhipeng Lu, 2015-10-16

Read in sam file with DG (duplex group) and XG (chiastic) tags,
filter out DGs with only 1 read and those with identical breaks, XG:i:2

"""

import sys, re, time
if len(sys.argv) < 3:
    print "samfilterdg.py"
    print "removes DGs with only single or duplicate reads and XG:i:2."
    print "Usage: python samfilterdg.py inputsam outputsam"
    sys.exit()

inputsam = sys.argv[1]
outputsam = sys.argv[2]
inputsamh = open(inputsam, 'r')
outputsamh = open(outputsam, 'w')

samheader = ''
numinputreads = 0
numoutputreads = 0
numoutputdg = 0
dgdict = {} #dgname: [dgreads]
for line in inputsamh: #construct a dictionary with all DGs
    if line[0] == "@":
        samheader += line
        continue
    numinputreads += 1
    record = line.strip('\n').split()
    if len(record) < 21: continue
    cigar = record[5]
    dgname = record[20]
    md = record[16].split(":")[-1]
    mdnum = len(re.findall('\d+[ATCG]', md))
    if record[19] == "XG:i:2" or len(dgname) <= 5 or mdnum > 1 : continue
    #remove reads with >1 mismatches, remove wrong DGs
    if not dgname in dgdict:
        dgdict[dgname] = [line]
    else: dgdict[dgname].append(line)
    if not numinputreads%10000:
        print time.strftime("%Y-%m-%d:%H:%M:%S"), "processed", numinputreads

#remove DG where all the reads have identical breaks in the CIGAR strings
#simply compare the N substrings for now
dgnamelist = dgdict.keys()
allreads = ''
for dgname in dgnamelist:
    dgreads = dgdict[dgname]
    breaklist = []
    for read in dgreads:
        record = read.strip('\n').split()
        cigar = record[5]
        cigarbits = tuple(re.findall('\d+[N]', cigar))
        breaklist.append(cigarbits)
    if len(list(set(breaklist))) == 1 :
        dgdict.pop(dgname)
        continue
    numoutputdg += 1
    numoutputreads += len(dgreads)
    dgout = ''.join(dgreads)
    allreads += dgout
    
outputsamh.write(samheader)
outputsamh.write(allreads)
print "\nNumber of input reads:", numinputreads
print "Number of filtered reads:", numoutputreads
print "Number of filtered duplex groups:", numoutputdg
inputsamh.close()
outputsamh.close()



"""
samheader = ''
numinputreads = 0
numoutputreads = 0
numoutputdg = 0
dgreadslist = []
outstring = ''

line = inputsamh.readline()
while line[0] == "@":
    samheader += line
    line = inputsamh.readline()
outputsamh.write(samheader)
record = line.strip('\n').split()
lastdgname = record[20]
dgreadslist.append(line)


for line in inputsamh:
    numinputreads += 1
    record = line.strip('\n').split()
    cigar = record[5]
    if len(record) <21: continue
    dgname = record[20]
    if dgname == lastdgname: dgreadslist.append(line)
    else:
        md = record[16].split(":")[-1]
        if "A" in md or "T" in md or "C" in md or "G" in md: continue
        breaklist = []
        for read in dgreadslist:
            record = read.strip('\n').split()
            cigar = record[5]
            cigarbits = tuple(re.findall('\d+[N]', cigar))
            breaklist.append(cigarbits)
        if len(list(set(breaklist))) > 1 :
            print "Passed filter:", dgname
            numoutputdg +=1
            numoutputreads += len(dgreadslist)
            dgout = ''.join(dgreadslist)
            outstring += dgout
        lastdgname = dgname
        dgreadslist = [line]
    if not numoutputdg%1000:
        outputsamh.write(outstring)
        outstring = ''
    
print "\nNumber of input reads:", numinputreads
print "Number of filtered reads:", numoutputreads
print "Number of filtered duplex groups:", numoutputdg
inputsamh.close()
outputsamh.close()
"""


