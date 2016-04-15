"""
alternativeconservation.py
This script examines the conservation/covariation signal from the list of
alternative structures. Then use Rchie to plot the conserved/covaried
alternative structures. Pay close attention to the DG name in each file!!!

Example command:
cd /Users/lu/Documents/chang/psoralen/higherorder
python ~/Documents/scripts/duplex/alternativeconservation.py \
~/Documents/chang/psoralen/examples/XIST/martinsmith/\
l15p2_NGmin1386_fix_COV/\
AMT_Stress_trim_nodup_starhsXIST_l15p2_NGmin1386_fix_Z2.326.bed \
AMT_Stress_trim_nodup_starhsXIST_l15p2_NGmin1386_fix.alt a
Result: 91 pairs out of 3462 pairs of XIST alternative structures (HEK)

cd /Users/lu/Documents/chang/psoralen/higherorder
python ~/Documents/scripts/duplex/alternativeconservation.py \
HeLa67HEK293_bed12_all_homo_sapiens_group_ge10sig.out \
AMT_Stress_trim_nodup_bc07_starhg38_geometricfilterfix_top50.alt a
Result: 7 pairs out of 448 pairs of alt structures in HEK293_1 top 50 mRNAs

cd /Users/lu/Documents/chang/psoralen/higherorder
python ~/Documents/scripts/duplex/alternativeconservation.py \
HeLa67HEK293_bed12_all_homo_sapiens_group_ge10sig.out \
AMT_HeLa6_trim_nodup_starhg38_geometricfilterfix_top50.alt a
Result: 31 pairs out of 711 pairs of alt structures in HeLa6 top 50 mRNAs

"""

import sys
if len(sys.argv) < 4:
    print "Usage: python alternativeconservation.py conserved alternative"
    print "       outlist"
    sys.exit()

conserved = open(sys.argv[1], 'r')
alternative = open(sys.argv[2], 'r')
outlist = open(sys.argv[3], 'w')

dgconserved = []
"""
for line in conserved: #if the input is in bed format
    record = line.split()
    dgname = record[3]
    dgconserved.append(dgname)
"""
for line in conserved: #if the input is in conservation format (custom)
    record = line.split()
    dgname = '_'.join(record[0].split('_')[0:2])
    dgconserved.append(dgname)

print dgconserved[0:3]


    
altconserved = 0
for line in alternative:
    record = line.split()
    altdg1 = '_'.join(record[3].split('_')[1:3])
    altdg2 = '_'.join(record[16].split('_')[1:3])
    #print altdg1, altdg2
    if altdg1 in dgconserved and altdg2 in dgconserved:
        altconserved += 1
        outlist.write(line)
        print line

print "Conserved alternative structure pairs:", altconserved

#Take the covariation information as well. 

    
"""
hsXIST 7710 7761 Group_1197 11 * 7710 7710 0,0,0 2 25,20 0,31
......((((.(.((((..(((((.&..))))))))).)))))...
hsXIST 7746 7885 Group_1223 3  * 7746 7746 0,0,0 2 20,20 0,119
(((((..((.((((((....&.....)))))).))))))).       0.75
"""


conserved.close()
alternative.close()
outlist.close()





