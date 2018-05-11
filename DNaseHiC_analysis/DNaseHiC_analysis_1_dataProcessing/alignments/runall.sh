#!/bin/bash -ex

# This script consists of two steps:
# STEP 1: bwa-mem alignment
# STEP 2: read filtering for primary and unmapped reads

STEP1=1
STEP2=1

DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"
BINDIR="<binariesDirName"

export PATH=$PATH:"${BINDIR}"/UCSC.utilities
assembly=mm10
if [ ! -e $assembly.sizes ]; then
	fetchChromSizes $assembly > $assembly.sizes
fi

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat libIDs)) # libIDs to array
IFS=$old_IFS

BWAIndexDir="$REFDIR}"/mm10pseudoSpretusBWAindex
assemblies=(black6 spretus)

for (( i=0; i<${#libIDs[@]}; i++ )); do
	libID=${libIDs[$i]}
	for (( j=0; j<${#assemblies[@]}; j++ )); do

		assembly=${assemblies[$j]}
		indices=$BWAIndexDir/$assembly.fa.gz

		if [[ "${STEP1}" -eq 1 ]]; then
			./step1_mapping.bwamem.sh $libID $indices $assembly $DATADIR
		fi 

		if [[ "${STEP2c}" -eq 1 ]]; then
			./step2_mappingSort.primaryAndUnmapped2sam.bwamem.sh $libID $assembly $DATADIR
		fi

	done
done

chmod 740 *.job
