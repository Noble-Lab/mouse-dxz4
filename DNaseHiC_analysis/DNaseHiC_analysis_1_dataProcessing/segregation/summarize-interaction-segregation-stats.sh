#!/bin/bash -x
#########################################################################
#
#code that generates the percentages of reads that:
# map to locations: =<200bp
#					 between 200bp 1kb
#					 >1kb
#					 inter
#
# 160315 Giancarlo Bonora: 
# 1. Made at least 6X faster 
# 2. added 1-20kb range
#
#########################################################################

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs.vSNPs)) # libIDs to array
IFS=$old_IFS

DATADIR=$(cat ../DATADIR)
CLEANEDDATADIR=${DATADIR}/$1

echo -e "LIBRARY\tend1\tend2\ttotal\tintra-close (<200bp)\tintra-mid (>=200bp & <1kbp)\tintra-20kb (>=1kbp & <20kbp)\tintra-far (>=20kbp)\tinter"

for i in ${libIDs[@]}; do
	thisLib=${CLEANEDDATADIR}/${i}/L1.gz
	for genome1 in ref alt both-ref; do
		for genome2 in ref alt both-ref; do
			if [[ ! "$genome1" > "$genome2" ]]; then
	zcat ${thisLib} | awk -v thisLib=${i} -v genome1=$genome1 -v genome2=$genome2 \
'$5==genome1 && $10==genome2 || $10==genome1 && $5==genome2 \
{NRgg++; if($3==$8) {x=$4-$9; \
if(x<200 && x>-200) {bpless200++}\
if((x>=200 && x<1000) || (x<=-200 && x>-1000)) {bp200v1kb++}\
if((x>=1000 && x<20000) || (x<=-1000 && x>-20000)) {bp1v20kb++}\
if(x<=-20000 || x>=20000) {bpmore20kb++}\
} else {bpinter++}}\
END{printf("%s\t%s\t%s\t%d (%.1f%%)\t\
%d (%.1f%%)\t\
%d (%.1f%%)\t\
%d (%.1f%%)\t\
%d (%.1f%%)\t\
%d (%.1f%%)\n",\
thisLib, genome1, genome2, NRgg, 100*NRgg/NR,\
bpless200, 100*bpless200/NRgg,\
bp200v1kb, 100*bp200v1kb/NRgg,\
bp1v20kb, 100*bp1v20kb/NRgg,\
bpmore20kb, 100*bpmore20kb/NRgg,\
bpinter, 100*bpinter/NRgg)}'
			fi
		done
	done
done
