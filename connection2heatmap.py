"""
connection2matrix.py
This program counts the connection 2xn data to a horizontal heatmap like the Hi-C data.
Note: when the matrix is big, plotting takes more time
example for hsXIST
python connection2matrix.py connectionfile matrixfile 50 19281

test:
python ~/Documents/scripts/paris/connection2heatmap.py AMT_Stress_trim_nodup_starhsXISTAligned_N.out.scatter AMT_Stress_trim_nodup_starhsXISTAligned_N.out_heatmap.pdf 100 19281
"""

import sys, math, numpy

if len(sys.argv) < 5: 
    print "Usage: python connection2heatmap.py connectionfile heatmapfile binsize rnasize"
    sys.exit()

connectionfile = sys.argv[1]
heatmapfile = sys.argv[2]
binsize = int(sys.argv[3])
rnasize = int(sys.argv[4])
connectionf = open(connectionfile, 'r')

r = math.sqrt(2)/2
numcol = rnasize*r*2/binsize+1
numrow = rnasize*r/binsize+1
matrix = numpy.zeros((numrow, numcol))
rotation = [[r, -r], [r, r]]


for line in connectionf:
    y, x = line.strip("\n").split()
    x, y = int(x), int(y)
    (x, y) = numpy.dot(rotation, [x, y])
    matrix[x/binsize, y/binsize] += 1

#rotmatrix = numpy.dot(rotation, matrix)
#print rotation

"""
a good resource:
http://www.labri.fr/perso/nrougier/teaching/matplotlib/
"""
import pylab, matplotlib.pyplot
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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



connectionf.close()
matrixf.close()
