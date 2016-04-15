import sys
import numpy as np
import matplotlib.pyplot as plt

xistcovfile = "/Users/lu/Documents/chang/psoralen/examples/XIST/martinsmith/\
l15p2_NGmin1386_fix_COV/AMT_Stress_trim_nodup_starhsXIST_l15p2_NGmin1386_fix.out"  

zscore = []
maskalncov = []
xistcov = open(xistcovfile, 'r')
for line in xistcov:
    if line[0:5] != "Group": continue
    record = line.split()
    zscore.append(float(record[8]))
    maskalncov.append(float(record[12]))
plt.scatter(zscore, maskalncov, alpha=0.5)
plt.xlim(-20, 2)
plt.savefig("zscore_maskalncov.pdf")
plt.show


"""
N = 50
x = np.random.rand(N)
y = np.random.rand(N)
colors = np.random.rand(N)
area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses

plt.scatter(x, y, s=area, c=colors, alpha=0.5)
plt.show()
"""
