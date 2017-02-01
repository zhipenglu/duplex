"""
converts duplex group data to bed12, bed12fixed or bed formats

example command for making 100nt intervals based on center of each arm
python ~/Documents/scripts/paris/dg2bed/dg2bed.py AMT_Stress_trim_nodup_starhsXIST_group_all_l15p2.support_filterRG AMT_Stress_trim_nodup_starhsXIST_group_all_l15p2.support_filterRG.bed bed12

Next steps:
1. use 
2. 

"""
import sys

if len(sys.argv) < 4:
    print "Usage: python dg2bed.py dgfile bedfile option\n\
           options as follows:\n\
           bed: bed6 format with 100nt intervals centered on each arm\n\
           bed12: regular bed12 format\n\
           bed12fixed: 40nt intervals each side\n"
    sys.exit()
    
dgfile = sys.argv[1]
bedfile = sys.argv[2]
option = sys.argv[3] #bed12 or bed12fixed or bed (each side separate)
                      #currently use 20nt each side for bed12fixed

#The problem with getting the proper region for phylogenetic analysis:
#If the two sides are too close: try extend more outside 
#The conservation/covariation methods do not restrict the base pairing between the two sides.
#

dgf = open(dgfile, 'r') # duplex group file
bedf = open(bedfile, 'w')
bedrecords = ''

#example duplex group:
#Group 564321 == position chr2(-):150469132-150469150|chr2(-):150487354-150487375, support 6/6.        

for line in dgf:
    group = line.split()
    if group and (group[0] == "Group"):
        bedinfo = group[4].replace("(", " ").replace(")", " ").replace(":", " ").replace("|", " ").split()
        chrom = bedinfo[0]
        bed1, bed2 = bedinfo[2].split("-")
        bed3, bed4 = bedinfo[5].strip(",").split("-")
        if int(bed3) < int(bed1):
            bed1, bed2, bed3, bed4 = bed3, bed4, bed1, bed2

        name = group[0] + "_" + group[1] + "_" + group[-1].strip(".")
        support = group[6].strip(',').split()[0] #this support is the number of reads
        strands = group[4].split("(")[1][0] + group[4].split("(")[2][0] 
        #"*" bedinfo[1] #strand info in the helix bed file makes visualization difficult. 
        itemRgb = "0,0,0"
        blockCount = "2"

        
        if option == "bed12": #for bed12:
            if strands[0] != strands[1]: continue
            strand = strands[0]
            chromStart = bed1
            chromEnd = bed4
            thickStart = chromStart
            thickEnd = chromStart
            blockSizes = str(int(bed2)-int(bed1)+1) + "," + str(int(bed4)-int(bed3)+1)
            blockStarts = "0," + str(int(bed3)-int(bed1))
            bedlist = [chrom, chromStart, chromEnd, name, support, \
                       strand, thickStart, thickEnd, itemRgb, \
                       blockCount, blockSizes, blockStarts]
            bedrecord = "\t".join(bedlist) + "\n"

        elif option == "bedpe": pass
            
            
        elif option == "bed12fixed": #for bed12fixed:
            if (int(bed3)+int(bed4))/2 - (int(bed1)+int(bed2))/2 < 20 :
                bed1new = str((int(bed1)+int(bed2)+int(bed3)+int(bed4))/4 - 20)
                bed2new = str((int(bed1)+int(bed2)+int(bed3)+int(bed4))/4)
                bed3new = str((int(bed1)+int(bed2)+int(bed3)+int(bed4))/4)
                bed4new = str((int(bed1)+int(bed2)+int(bed3)+int(bed4))/4 + 20)
                bed1, bed2, bed3, bed4 = bed1new, bed2new, bed3new, bed4new
            else:
                bed1new = str((int(bed1)+int(bed2))/2 - 10)
                bed2new = str((int(bed1)+int(bed2))/2 + 10)
                bed3new = str((int(bed3)+int(bed4))/2 - 10)
                bed4new = str((int(bed3)+int(bed4))/2 + 10)
                bed1, bed2, bed3, bed4 = bed1new, bed2new, bed3new, bed4new
            chromStart = bed1
            chromEnd = bed4
            thickStart = chromStart
            thickEnd = chromStart
            blockSizes = str(int(bed2)-int(bed1)+1) + "," + str(int(bed4)-int(bed3)+1)
            blockStarts = "0," + str(int(bed3)-int(bed1))
            bedlist = [chrom, chromStart, chromEnd, name, support, \
                       strand, thickStart, thickEnd, itemRgb, \
                       blockCount, blockSizes, blockStarts]
            bedrecord = "\t".join(bedlist) + "\n"

        
        elif option == "bed": #for simple bed with 100nt interval:
            leftcenter = (int(bed1)+int(bed2))/2
            rightcenter = (int(bed3)+int(bed4))/2
            bedleft = []
            bedright = []
            if int(support) >= 25:
                bedleft = [chrom, str(leftcenter-50), str(leftcenter+50), support]
                bedright = [chrom, str(rightcenter-50), str(rightcenter+50), support]
            bedrecord = "\t".join(bedleft) + "\n" + "\t".join(bedright) + "\n"


        bedrecords += bedrecord
        
bedf.write(bedrecords)

dgf.close()
bedf.close()
