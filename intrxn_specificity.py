#snorna_specificity.py
#this script plots the specificity of RNA:RNA basepairing
#basic steps:
#read the interactions from
#AMT_Stress_trim_nodup_norm_starRfamhumanrna_geometric.bed or
#
#then plot the coverage for one RNA of the interaction

import sys, re
import matplotlib.pyplot as plt
from pyliftover import LiftOver
if len(sys.argv) < 1:
    print "Usage: python snorna_specificity.py"
    print "change the parameters in the script"
    sys.exit()


#parameters [intrxnfile, RNAtoplot, RNAsize, asp_ratio, outpdf, xlab, ylab]
SNORD14_on_18S_HEK = ["SNORD14", "hs18S", 1869, "SNORD14_on_18S_HEK.pdf", \
                      "human 18S rRNA", "SNORD14:hs18S \n coverage (HEK)"]
SNORD95_on_28S_HEK = ["SNORD95", "hs28S", 5070, "SNORD95_on_28S_HEK.pdf", \
                      "human 28S rRNA", "SNORD95:hs28S \n coverage (HEK)"]
SNORD95_on_28S_mES = ["mm28S", "SNORD95", 5070, "SNORD95_on_28S_mES.pdf", \
                      "mouse 28S rRNA (human coordinates)",\
                      "SNORD95:mm28S \n coverage (mES)"]
snou13_on_18S_HEK =  ["snoU13", "hs18S", 1869, "SNOU13_on_18S_HEK.pdf", \
                      "human 18S rRNA", "SNOU13:18S \n coverage (HEK)"]
snoU13_on_18S_mES =  ["snoU13", "mm18S", 1869, "SNOU13_on_18S.pdf", \
                      "mouse 18S rRNA", "SNOU13:18S \n coverage (mES)"]
SNORD10_on_U6_HEK = ["SNORD10", "U6", 107, "SNORD10_on_U6_HEK.pdf", \
                     "human U6 snRNA", "SNORD10:U6 \n coverage (HEK)"]
SNORD10_on_U6_mES = ["SNORD10", "U6", 108, "SNORD10_on_U6_mES.pdf", \
                     "mouse U6 snRNA", "SNORD10:U6 \n coverage (mES)"]
U3_on_45S_HeLa =  ["U3", "hs18S", 13357, "U3_on_45S_HeLa.pdf", \
                   "human 45S rRNA", "U3:18S \n coverage (HeLa)"]
U3_on_45S_HEK =  ["U3", "hs45S", 13357, "U3_on_45S_HEK.pdf", \
                  "human 45S rRNA", "U3:45S \n coverage (HEK)"]
U3_on_45S_mES =  ["U3", "mm45S", 13357, "U3_on_45S_mES.pdf", \
                  "mouse 45S rRNA (human coordinates)", "U3:45S \n coverage (mES)"]

U8_on_45S_HEK =  ["U8", "hs45S", 13357, "U8_on_45S_HEK.pdf", \
                  "human 45S rRNA", "U8:45S \n coverage (HEK)"]
U8_on_45S_mES =  ["U8", "mm45S", 13357, "U8_on_45S_mES.pdf", \
                  "mouse 45S rRNA (human coordinates)", "U8:45S \n coverage (mES)"]
hs28S_on_U8_HEK =  ["hs28S", "_U8", 135, "hs28S_on_U8_HEK.pdf", \
                    "human U8 snoRNA", "U8:28S \n coverage (HEK)"]
mm28S_on_U8_mES =  ["hs28S", "_U8", 135, "mm28S_on_U8_mES.pdf", \
                    "mouse U8 snoRNA", "U8:28S \n coverage (mES)"]
RN7SK_on_28S_HEK = ["7SK", "hs28S", 5070, "7SK_on_28S_HEK.pdf", \
                    "human 28S rRNA", "7SK:28S \n coverage (HEK)"]
RN7SK_on_28S_mES = ["7SK", "mm28S", 5070, "7SK_on_28S_mES.pdf", \
                    "mouse 28S rRNA (human coordinates)", "7SK:28S \n coverage (mES)"]
RN7SK_on_18S_HEK = ["7SK", "hs18S", 1869, "7SK_on_18S_HEK.pdf", \
                    "human 18S rRNA", "7SK:18S \n coverage (HEK)"]
RN7SK_on_18S_mES = ["7SK", "mm18S", 1870, "7SK_on_18S_mES.pdf", \
                    "mouse 18S rRNA", "7SK:18S \n coverage (mES)"]
snoU83B_on_18S_HEK = ["snoU83B", "hs18S", 1869, "snoU83B_on_18S_HEK.pdf", \
                      "human 18S rRNA", "snoU83B:18S \n coverage (HEK)"]
snoU83B_18S_mES =  ["snoU83B", "mm18S", 1870, "snoU83B_on_18S_mES.pdf", \
                    "mouse 18S rRNA", "snoU83B:18S \n coverage (mES)"]
snoU83B_28S_HEK = ["snoU83B", "hs28S", 5070, "snoU83B_on_28S_HEK.pdf",\
                   "human 28S rRNA", "snoU83B:28S \n coverage (HEK)"]
snoU83B_28S_mES =  ["snoU83B", "mm28S", 5070, "snoU83B_on_28S_mES.pdf", \
                    "mouse 28S rRNA (human coordinates)", \
                    "snoU83B:28S \n coverage (mES)"]
U1_on_XIST_HEK =  ["_U1", "XIST", 19285, "U1_on_XIST_HEK.pdf", \
                   "human XIST", "U1:XIST \n coverage (HEK)"]
U1_on_MALAT1_HEK =  ["_U1", "MALAT1", 8708, "U1_on_MALAT1_HEK.pdf",
                     "human MALAT1", "U1:MALAT1 \n coverage (HEK)"]
U1_on_Malat1_mES =  ["Malat1", "_U1", 8708, "U1_on_Malat1_mES.pdf", \
                     "mouse Malat1 (human coordinates)", "U1:Malat1 \n coverage (mES)"]
MALAT1_on_U1_HEK =  ["MALAT1", "U1", 165, "MALAT1_on_U1_HEK.pdf", \
                     "human U1", "U1:MALAT1 \n coverage (HEK)"]
MALAT1_on_U1_mES =  ["Malat1", "U1", 165, "MALAT1_on_U1_mES.pdf", \
                     "mouse U1 (human coordinates)", "U1:Malat1 \n coverage (mES)"]
SNORD16_on_U6_HEK = ["SNORD16", "_U6", 107, "SNORD16_on_U6_HEK.pdf", \
                     "human U6", "SNORD16:U6 \n coverage (HEK)"]
SNORD16_on_U6_mES = ["SNORD16", "_U6", 107, "SNORD16_on_U6_mES.pdf", \
                     "mouse U6", "SNORD16:U6 \n coverage (mES)"]
U6_on_SNORD16_HEK = ["_U6", "SNORD16", 99, "U6_on_SNORD16_HEK.pdf", \
                     "human SNORD16", "SNORD16:U6 \n coverage (HEK)"]
U6_on_SNORD16_mES = ["_U6", "SNORD16", 99, "U6_on_SNORD16_mES.pdf", \
                     "mouse SNORD16", "SNORD16:U6 \n coverage (mES)"]
SNORD101_on_28S_HEK = ["SNORD101", "28S", 5070, "SNORD101_on_28S_HEK.pdf", \
                       "human 28S rRNA", "SNORD101:28S \n coverage (HEK)"]
SNORD101_on_28S_mES = ["SNORD101", "mm28S", 5070, "SNORD101_on_28S_mES.pdf", \
                       "mouse 28S rRNA (human coordinates)", "SNORD101:28S \n coverage (mES)"]
SNORA64_on_28S_HEK = ["SNORA64", "hs28S", 5070, "SNORA64_on_28S_HEK.pdf", \
                      "human 28S rRNA", "SNORA64:28S \n coverage (HEK)"]
SNORA73_on_28S_HEK = ["SNORA73", "hs28S", 5070, "SNORA73_on_28S_HEK.pdf", \
                      "human 28S rRNA", "SNORA73:28S \n coverage (HEK)"]
SNORA73_on_28S_mES = ["SNORA73", "mm28S", 5070, "SNORA73_on_28S_mES.pdf", \
                      "mouse 28S rRNA", "SNORA73:28S \n coverage (mES)"]
SNORA73_on_18S_HEK = ["SNORA73", "hs18S", 1870, "SNORA73_on_18S_HEK.pdf", \
                      "human 18S rRNA", "SNORA73:18S \n coverage (HEK)"]
SNORA73_on_18S_mES = ["SNORA73", "mm18S", 1870, "SNORA73_on_18S_mES.pdf", \
                      "mouse 18S rRNA", "SNORA73:18S \n coverage (mES)"]
SNORD22_on_28S_HEK = ["SNORD22", "hs28S", 5070, "SNORD22_on_28S_HEK.pdf", \
                      "human 28S rRNA", "SNORD22:28S \n coverage (HEK)"]
SNORD22_on_28S_mES = ["SNORD22", "mm28S", 5070, "SNORD22_on_28S_mES.pdf", \
                      "mouse 28S rRNA", "SNORD22:28S \n coverage (mES)"]



HEK_intrxnfile = "AMT_Stress_trim_nodup_norm_starRfamhumanrnaMrna_geometricfiltered"
mES_intrxnfile = "AMT_Mettl3_trim_nodup_wt_starRfammousernaMrna_geometricfiltered"
current = U3_on_45S_mES
intrxnfile = open(mES_intrxnfile, 'r')
partner = current [0]
RNAtoplot = current[1]
size = current[2]
outpdf = current[3]
xlab, ylab = current[4], current[5]


dist = [0 for i in range(0, size)] # initialize a list to store coverage 
#Use this file as input:
#AMT_Stress_trim_nodup_norm_starRfamhumanrnaMrna_geometricfiltered
#example line: readname 7SK|+:64-86<=>PLCG2|+:7967-8004
rightpair = 0
for line in intrxnfile:
    if line[0] == "G": continue
    intrxn = line.strip('\n').split()[1].replace("<=>", " ").split()
    plotrna, interval = '', ''
    if RNAtoplot in intrxn[0] and partner in intrxn[1]:
        interval = intrxn[0].split(":")[1].split("-")
    elif RNAtoplot in intrxn[1] and partner in intrxn[0]:
        interval = intrxn[1].split(":")[1].split("-")
    if len(interval) == 2:
        for i in range(int(interval[0]), int(interval[1])):
            dist[i] += 1
print "RNA size:", len(dist)


#Use the following part to liftover mouse coordinates to human 
liftfiles = {"mm28S": "/Users/lu/Documents/chang/rrna/liftover/mmtohs28S.liftoverchain", \
"mm45S": "/Users/lu/Documents/chang/rrna/liftover/mmtohs45S.liftoverchain", \
"Malat1": "/Users/lu/Documents/chang/psoralen/examples/MALAT1/mmtohg_Malat1.liftoverchain"}
if RNAtoplot in liftfiles:
    newdist = [0 for i in range(0, size)]
    lo = LiftOver(liftfiles[RNAtoplot])
    for i in range(0, size):
        lifted = lo.convert_coordinate(RNAtoplot, i, '+')
        if lifted: newdist[lifted[0][1]] += dist[i]
    dist = newdist



figure = plt.figure(figsize=(8,2))
axes = plt.Axes(figure, [.3,.3,.6,.6])
figure.add_axes(axes)
plt.bar(range(0, size), dist, color='k')
axes.spines['top'].set_visible(False)
axes.spines['right'].set_visible(False)
axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
plt.xlim(0, size)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.savefig(outpdf)
plt.show()



