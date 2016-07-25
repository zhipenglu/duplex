#!/bin/bash
# Requires Jim Kent's mafsInRegions binary from UCSC tools
# Requires bedtools 
# Requires Vienna RNA version >1.8.5
# Requires SISSIz version >2


###################################################################
###################################################################
#change the following to fit this computer: 
#remove the RNAalifold path
#change RNAalifold options from -color to --color, and from -aln to --aln
#when Convert to Clustal format, add the $1=="s" condition to the awk command
#make sure that the chromosome names are correct: hg or homo_sapiens, depending on the maf chr names. 
#example: sed 's/chr/homo_sapiens./g' PFN1_AMT_HeLa6_trim_nodup_starhg38_l15p2_geometric_NGmin_fix.bed > a.bed
#./PFN1_covariation.sh PFN1_AMT_HeLa6_trim_nodup_starhg38_l15p2_geometric_NGmin_fix.bed | grep ^Group > MALAT1_AMT_Stress_trim_nodup_norm_starhg38_l15p2_geometric_filteredfix.out

#only need to change the maf file in the get alignments section. 
###################################################################
###################################################################







########################################
# Format alignment and convert mafft-ginsi fasta to MAF
########################################
#awk 'OFS="\t" { if ( $1 == "s" ) print $1,$2,$3,$4,$5,($6+1),$7 ; else print }' hsXIST_multiz100_exons_eutherian58.maf |\
#	sed 's/hsXIST/hg\.hsXIST/g' > hsXIST_multiz100_exons_eutherian58_f.maf
#( echo "##maf version=1.0" 
#  echo "a score=1"
#  awk '{ if (NR == 1 ) printf $1"\t"; else if ( match( $1 , ">" )) printf "\n"$1"\t" ; else printf $1 }'  hsXIST_mafft.fasta |\
#	 cut -c 2- | join -1 1 -2 2  -o 2.1,1.1,2.3,2.4,2.5,2.6,1.2 - <( tail -n +3  hsXIST_multiz100_exons_eutherian58_f.maf) | awk ' OFS="\t"{print}' 
#) > hsXIST_multiz100_exons_eutherian58_mafft.maf
#
########################################
# process BEDs
########################################
# N.B. "*" as strand inverses the order of the sequences (2 is 5'/upstream, of 1)
join -1 4  <( cat $1 | bed12ToBed6 -n | sort -k 4 ) <( cut -f 4-5 $1 | sort -k 4,4 ) |\
 awk 'OFS="\t" {print $2,$3,$4,$1"_"$5,$7,$6}'> ${1%*.bed}_6.bed








########################################
# get alignments using maf_extract_ranges_indexed.py which is much faster than mafsInRegion
########################################
if [[ ! -d ./maf_chunks ]]; then mkdir ./maf_chunks ; fi 
cat ${1%*.bed}_6.bed | while read line ; do 
echo $line | python ~/Documents/analysis/bxpython/scripts/maf_extract_ranges_indexed.py -c $2 > maf_chunks/$(grep -o 'Group\S*' <<< $line).maf
done

cd ./maf_chunks
for file1 in *_1.maf; do
file2=`echo $file1 |sed 's/_1.maf/_2.maf/g'`
line1=`cat $file1 | wc -l`
line2=`cat $file2 | wc -l`
count1=`awk '($1=="s")&&($4>100)' $file1 |wc -l`
count2=`awk '($1=="s")&&($4>100)' $file2 |wc -l`
if [ $((count1+count2)) -gt 0 -o $line1 -lt 8 -o $line2 -lt 8 ]; then rm $file1 $file2; fi
done
cd ../

# Convert to Clustal format
for file in ./maf_chunks/*maf 
do 
 (echo -e "CLUSTAL\n"; tail -n +3 $file | awk '$1 == "s" { printf("%-25s %s\n",$2,$7) }' ) > ${file%*maf}aln
done
 
########################################
# make a mock-hairpin for RNAalifold
########################################
if [[ ! -d ./mock_hairpins ]]; then mkdir ./mock_hairpins ; fi 
	cut -f 4 $1 | while read line ; do join ./maf_chunks/${line}_2.aln ./maf_chunks/${line}_1.aln |\
  awk '{ if (NF > 1 ) printf("%-25s %s\n",$1,$2"----------"$3); else print }' | awk '{ if ( $2 !~ /^[-]+$/ ) print }' > ./mock_hairpins/${line}_hairpin.aln; done

########################################
# Run RNAalifold and SISSIz on mock hairpins
# Also get RNAcofold output for constrained and unconstrained duplexes
########################################
# print header
echo -e -n 	"#File\tsize\tAln%ID\tSim%ID\tSim%ID-stdev\tAln_MFE\tSim_MFE_mean\tSim_MFE_stdev\tZ-score\t"
echo -e 	"Aln_MFE\tAln_COV\tMaskAln_MFE\tMaskAln_COV\tcofold_naive\tcofold_constr\tnaive_2D\tconstr_2D\tConsensus_Seq\tConsensus_2D"
if [[ ! -d ./figs ]]; then mkdir ./figs ; fi 
cd mock_hairpins
for file in *hairpin.aln
do 
	##cofold 
	##	echo `pwd`
	SEQ2=$( grep hs45S ../maf_chunks/${file%_hairpin.aln}_2.aln | awk '{print $2}' | sed s/-//g )
	echo $SEQ2
	SEQ1=$( grep hs45S ../maf_chunks/${file%_hairpin.aln}_1.aln | awk '{print $2}' | sed s/-//g ) 
	STR2=$( echo $SEQ2 | sed 's/[ACGUTacgut-]/</g') 
	STR1=$( echo $SEQ1 | sed 's/[ACGUTacgut-]/>/g')
	echo $SEQ2"&"$SEQ1
	COFOLD1=$( echo -e $SEQ2"&"$SEQ1 | RNAcofold  | tail -n +1 )
	mv -f ./rna.ps ../figs/${file%*.aln}_COFOLD_naive.ps
	COFOLD2=$( echo -e $SEQ2"&"$SEQ1"\n"$STR2"&"$STR1 | RNAcofold -C  | tail -n +1 )
	mv -f ./rna.ps ../figs/${file%*.aln}_COFOLD_constrained.ps
	MFE1=$( echo $COFOLD1 | cut -d " " -f 3- | sed -e 's/[\(\)]//g' -e 's/ +//g' )
	MFE2=$( echo $COFOLD2 | cut -d " " -f 3- | sed -e 's/[\(\)]//g' -e 's/ +//g' )

	##alifold - N.B. newer alifold parameters are --color --aln (double dash)
	OUT=$( RNAalifold --color --aln -r ${file}  2>/dev/null ) 
	mv -f ./alirna.ps ../figs/${file%*.aln}_RNA2D.ps
	mv -f ./aln.ps ../figs/${file%*.aln}_ALN.ps
	##Constrained
	MASK2=$( grep homo_sapiens ../maf_chunks/${file%_hairpin.aln}_2.aln |  awk '{print $2}' | sed 's/[ACTGUactgu-]/</g' )
	MASK1=$( grep homo_sapiens ../maf_chunks/${file%_hairpin.aln}_1.aln |  awk '{print $2}' | sed 's/[ACTGUactgu-]/>/g' )
	echo $MASK1
	OUTC=$( echo ${MASK2}"----------"${MASK1} | RNAalifold --color --aln -C -r ${file}  2>/dev/null ) 
	#echo ${MASK2}"----------"${MASK1}
	#echo $file

	mv -f ./alirna.ps ../figs/${file%*.aln}_RNA2D_MASK.ps
	mv -f ./aln.ps ../figs/${file%*.aln}_ALN_MASK.ps

	echo -e "[ NOTE ] Running 100 simulations with SISSIz on "$file  >&2
	SIM=$( SISSIz -j -n 100 $file | cut -f 2,4- ) 
	echo -n $SIM" " 
	echo -n $OUT | cut -d " " -f 3- | sed -e s/\(// -e s/\)// -e s/=// -e s/+// | awk '{printf $2" "$3" "}' 
	echo -n $OUTC | cut -d " " -f 3- | sed -e s/\(// -e s/\)// -e s/=// -e s/+// | awk '{printf $2" "$3" "}' 
	echo -n $MFE1" "$MFE2" "`echo $COFOLD1 | cut -d " " -f 2`" "`echo $COFOLD2 | cut -d " " -f 2`" "
	echo -n $OUT | awk '{printf $1" "$2" "}'
	echo -n $OUTC | awk '{printf $2" \n"}'
done | awk 'OFS="\t" {print}'



# ########################################
# # OPTION: Run RNAalifold on mock haipins with inter-sequence constraints
# for file in *hairpin.aln
# do
# 	S2=$( grep hg ./maf_chunks/${file%*_hairpin.aln}_2.aln | awk '{print $2}' | sed 's/[ATCG-]/</g' )
# 	S1=$( grep hg ./maf_chunks/${file%*_hairpin.aln}_1.aln | awk '{print $2}' | sed 's/[ATCG-]/>/g' )
# 	echo $file 
# 	echo $S2".........."$S1 | RNAalifold --color --aln --mis -r -C $file 2>/dev/null | awk 'OFS="\t"{ if (NR==1) seq=$1; else if (NR == 2) {struct = $1; mfe = $5 ; cov = $7} } END {print mfe,cov,seq,struct}' 
# done
 
# ########################################
# # simulate alignments (shuffle each block independently) 
# # N.B. the simluation might work better if we shuffle the entire "hairpin" alignment
# ########################################

# cut -f 4 $1 | while read line 
# do 
#  #To do : add a while loop for significance
#  echo $line
# 	(echo -e "CLUSTAL   ${line}\n";
# 	join <(SISSIz -s -f 100 --maf ./maf_chunks/${line}_2.maf |\
# 		cut -d " " -f 2,7) <(SISSIz -s -f 100 --maf ./maf_chunks/${line}_1.maf |\ 
# 		cut -d " " -f 2,7) |\
# 		awk '{ if ( NF > 1 ) printf "%-25s%s\n",$1,$2"----------"$3 ; else print }' |\
# 		awk '{ if ( $2 !~  /^-*$/) print }' ) |\
# 		RNAalifold -r 
#  # To do : collect RNAalifold scores for Z-score calculation		
# done

