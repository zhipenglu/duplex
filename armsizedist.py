"""
this script calculates the size distribution of the size of each arm in both
the normal gapped read file: Aligned_prim_N.sam (two longest M each record)
and the chiascti file: Chimeric.out.sam (one longest M each record)

command for creating the list:
python ~/Documents/scripts/paris/armsizedist.py \
AMT_HeLa6_trim_nodup_starhg38Aligned_prim_N.out_sorted.sam \
AMT_HeLa6_trim_nodup_starhg38Chimeric.out_sorted.sam AMT_HeLa6_armsize.txt


"""

import sys, re, matplotlib
import matplotlib.pyplot as plt

if len(sys.argv)< 4:
    print "Usage: python armsizedist.py aligned.sam chimeric.sam outfile"
    #sys.exit()

#alignedsamfile = sys.argv[1]
#chimericsamfile = sys.argv[2]
#outfile = sys.argv[3]



"""
#quote out this section once the size list is created well.
#get armsize in alignedsamfile
sizelist = []
alignedsamfileh = open(alignedsamfile, 'r')
for line in alignedsamfileh:
    if line[0] == "@": continue
    sam = line.split()
    cigar = sam[5]
    cigars = re.findall('[0-9]+[M]', cigar)
    matches = [int(match.strip("M")) for match in cigars]
    sizelist += sorted(matches)[-2:]
alignedsamfileh.close()
#get armsize in chimericsamfile
chimericsamfileh = open(chimericsamfile, 'r')
for line in chimericsamfileh:
    if line[0] == "@": continue
    sam = line.split()
    cigar = sam[5]
    cigars = re.findall('[0-9]+[M]', cigar)
    matches = [int(match.strip("M")) for match in cigars]
    sizelist.append(sorted(matches)[-1])
chimericsamfileh.close()
#save the size list as space delimited numbers in a single line in a file. 
sizestring = ' '.join([str(size) for size in sizelist])
outfileh = open(outfile, 'w')
outfileh.write(sizestring)
outfileh.close()
"""


#plot the dg size distribution
outfileh = open("HEK293pe.dgsize", 'r')
sizelist = [int(i) for i in outfileh.readline().split()]
outfileh.close()
print len(sizelist)
fig, ax = plt.subplots()
n, bins, patches = plt.hist(sizelist, [x*10+0.5 for x in range(0,1000)], histtype='stepfilled', cumulative=True)#normed=1, 
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
plt.xlim(0, 10000)
ax.set_xlabel("DG size", fontsize=15)
ax.set_ylabel("frequency", fontsize=15)
plt.savefig("HEK293pe.dgsize.pdf")
plt.show()





"""
#plot the arm size distribution
outfileh = open("AMT_Stress_trim_nodup_starhsXIST_l15p2_filtered_support10.armsize", 'r')
sizelist = [int(i) for i in outfileh.readline().split()]
print len(sizelist)
fig, ax = plt.subplots()
n, bins, patches = plt.hist(sizelist, [x+0.5 for x in range(0,80)], histtype='stepfilled')#normed=1, 
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
plt.xlim(1, 80)
ax.set_xlabel("arm length (nt, gapped and chiastic, support >= 10)", fontsize=15)
ax.set_ylabel("frequency", fontsize=15)
plt.savefig("AMT_Stress_trim_nodup_starhsXIST_l15p2_filtered_support10.armsize.pdf")
plt.show()
"""



"""

plt.figure(facecolor="white")
fig = plt.pcolor(matrix, cmap="Blues", vmin=0, vmax=200)
plt.axis()
plt.axes().set_aspect('equal')
locs = range(0, int(rnasize/(binsize*r*2)), int(5000/(binsize*r*2)))
labels = range(0, rnasize, 5000)
plt.yticks(locs, labels)
plt.xlim(0, rnasize*r*2/float(binsize)+1)
plt.ylim(0, rnasize*r/float(binsize)+1)
plt.axes().spines['top'].set_visible(False)
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['bottom'].set_visible(False)
plt.axes().yaxis.set_ticks_position('left')
plt.axes().tick_params(top="off")
plt.axes().tick_params(bottom="off")
plt.axes().set_xticks([])

cb = plt.colorbar(mappable=fig,
                  orientation="horizontal",
                  shrink=0.3,
                  aspect=20,
                  ticks=[0,50,100,150,200]
                  )
cb.ax.set_title('PARIS coverage')
cb.ax.set_xticklabels(['0', '50', '100', '150', '>200'])
plt.savefig(heatmapfile)
plt.show()

"""


