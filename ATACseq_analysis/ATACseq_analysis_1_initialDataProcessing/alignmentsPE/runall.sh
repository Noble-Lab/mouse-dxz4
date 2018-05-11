#!/bin/bash -ex
#################################################
#
# 20170105
# Giancarlo
# BWA PE mapping to each assembly
#
#################################################

set -o nounset
set -o pipefail
set -o errexit 

DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"
BINDIR="<binariesDirName"

export PATH=$PATH:"${BINDIR}"/UCSC.utilities
assembly=mm10
if [[ ! -e $assembly.size ]]; then
	fetchChromSizes $assembly > $assembly.sizes
fi

MISMATCH=3

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs)) # libIDs to array
IFS=$old_IFS

#assemblies=(hg18 hg19)
#assemblies=(black6 spretus spretus.Fan.combined)
#assemblies=(black6 spretus.Fan.combined)

BWAIndexDir=$HOME/refdata/mm10pseudoSpretusBWAindex
assemblies=(black6 spretus)

STEP1=1
STEP2=1

for (( i=0; i<${#libIDs[@]}; i++ )); do
  libID=${libIDs[$i]}
    for (( j=0; j<${#assemblies[@]}; j++ )); do
      assembly=${assemblies[$j]}
      indices=$BWAIndexDir/$assembly.fa.gz

      if [[ "${STEP1}" -eq 1 ]]; then
        ./step1_mapping.bwamem.PE.sh $libID $indices $assembly $DATADIR
	  fi 
   
      if [[ "${STEP2}" -eq 1 ]]; then
	    ./step2_mappingSort.primaryAndUnmapped2sam.bwamem.PE.sh $libID $assembly $DATADIR
      fi
  done
done

chmod 740 *.job