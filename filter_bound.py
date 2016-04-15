#input file: sam files with or without headlines
#input range, no chromosome information needed
#Only "include" alignments with both ends within the range
#OR "exclude" those with both ends within the ranges.

#Note only "MINDS" are considered from the CIGAR string
#Author: Zhipeng Lu, zhipengluchina@gmail.com
#2015-10-02

"""
take a random sample from this file:
cd ~/Documents/chang/psoralen/HEK1/HEK1_hs45S/
python ~/Documents/scripts/duplex/randlines.py \
AMT_Stress_trim_nodup_bc07_starhs45SAligned_N.out_sorted.sam \
10000 hs45S_10000.sam

test command on the 45S rRNA:
cd ~/Documents/chang/psoralen/HEK1/HEK1_hs45S/
python ~/Documents/scripts/filter_bound.py hs45S_10000_sorted.sam precursor \
0,13357 3654,6600,7924 5523,6757,12994

python /home/zhipeng/bin/filter_bound.py \
AMT_Stress_trim_nodup_bc07_starhs45SAligned_N.out_sorted.sam precursor \
0,13357 3654,6600,7924 5523,6757,12994

python /home/zhipeng/bin/filter_bound.py \
AMT_Mettl3_trim_nodup_bc02_starmm45SAligned_N.out_sorted.sam precursor \
0,13400 4008,6878,8123 5877,7034,12849

example read:
#M01339    256     chr21   7814836 0       10M1I10M403206N39M93N85M18S
"""

import sys, re
if len(sys.argv) < 6:
    print "Usage: python filter_bound.py input.sam [include|exclude|precursor] \
region_bounds left_bounds right_bounds\n\
       left_bounds and right_bounds are comma-separated lists\n\
       Only 'include' alignments with both ends within the range(s) OR\n\
       'exclude' alignments with both ends within the range(s),\n\
       'precursor' means removing reads mapped entirely in the mature regions.\n\
       left/right_bounds indicate the mature regions of the RNA.\n\
       note: only [MINDS] operations from the CIGAR string are considered\n"
    sys.exit()

inputsam = sys.argv[1]
flag = sys.argv[2]
region_bounds = [int(i) for i in sys.argv[3].split(',')]
left_bounds  =  [int(i) for i in sys.argv[4].split(',')]
right_bounds =  [int(i) for i in sys.argv[5].split(',')]
removed_regions = []
left = left_bounds + region_bounds[1:]
right = region_bounds[0:1] + right_bounds
print right, left

if flag == "include":
    outsam = sys.argv[1].strip(".sam") + "_include.sam"
elif flag == "exclude":
    outsam = sys.argv[1].strip(".sam") + "_exclude.sam"
elif flag == "precursor":
    outsam = sys.argv[1].strip(".sam") + "_precursor.sam"
else: print "please use 'include' or 'exclude' as the flag"


inputf = open(inputsam, 'r')
outf = open(outsam, 'w')
for line in inputf:
    if line[0] == "@":
        outf.write(line)
    else:
        record = line.split()
        start = int(record[3])
        cigars = re.findall('\d+[MIDNSPH=X]', record[5])
        readrange = 0
        for cigar in cigars:
            if cigar[-1] in ["M", "D", "N"]: readrange += int(cigar[:-1])
            elif cigar[-1] in ["I"]: readrange -= int(cigar[:-1])
            #do nothing if the cigar operation is "S" or the others [HP=X].
        end = start + readrange
        
        if flag == "include":
            for i in range(len(left_bounds)):
                if (start >= left_bounds[i]) and (end <= right_bounds[i]):
                    outf.write(line)
        elif flag == "exclude":
            if end < left_bounds[0] or start > right_bounds[-1]:
                outf.write(line)
            else:
                for i in range(len(left_bounds)-1):
                    if (start >= right_bounds[i]) and (end <= left_bounds[i+1]):
                        outf.write(line)
        elif flag == "precursor":
            for i in range(len(left)):
                if start > right[i] and start < left[i] or \
                   end > right[i] and end < left[i]:
                    outf.write(line)
        else: print "please use 'include' or 'exclude' as the flag"
        
inputf.close()
outf.close()





