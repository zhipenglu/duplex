#rfam2rna.py
#this script retrieves the sequences for ncRNAs from a list of rfam ids.
#example command:
#python ~/Documents/scripts/rfam2rna.py Rfam.seed Rfam.humanrna Rfam_humanrna.fa


import sys

if len(sys.argv) < 4:
    print "Usage: python rfam2rna.py Rfam.seed Rfam.humanrna Rfam_humanrna.fa"
    sys.exit()

rfamseed = open(sys.argv[1], 'r')
rfamlist = open(sys.argv[2], 'r')
rfamrnafa = open(sys.argv[3], 'w')



#     format of the seed file
#     #=GF AC   RF00001
#     #=GF ID   5S_rRNA
#     X01556.1/3-118
rfamdict = {} #store the rfam seed file in a dictionary.
AC, ID = '', ''
familycount = 0
for line in rfamseed:
    record = line.strip('\n').split()
    if len(record) < 2: continue
    if record[1] == "AC":
        AC = record[2]
        familycount += 1 
    if record[1] == "ID": ID = record[2]
    if line[0] != "#":
        seqname = record[0].split('/')[0]
        seq = record[1].translate(None, '-').replace("U", "T")
        rfamdict[(AC, seqname)] = [ID, seq]
print "Number of families:", familycount
      
            

#read the rfam family and rna id information and extract sequence
outfa = ''
for line in rfamlist:
    record = line.strip('\n').split()
    AC, ID = record[0], record[1]
    seqname, seq = tuple(rfamdict[(AC, ID)])
    fa = ">" + AC + "_" + seqname + "\n" + seq + "\n"
    outfa += fa
rfamrnafa.write(outfa)

rfamlist.close()
rfamrnafa.close()
