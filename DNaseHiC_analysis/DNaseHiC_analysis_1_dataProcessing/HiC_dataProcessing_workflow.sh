##################################################
##################################################
##################################################
# 
# Map all Hi-C libraries and generate contact maps
#
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="<projDirName>"
mkdir -p "${PROJDIR}"
DATADIR="<dataDirName>"





##################################################
##################################################
##################################################
# Read alignment

mkdir "${PROJDIR}"/alignments
cd "${PROJDIR}"/alignments

##################################################
# Runs BWA
# ***       2 STEPS:
# ***       1) STEP 1: INITIAL BWA-MEM ALIGNMENT AND 
# ***       2) STEP 2: READ FILTERING FOR PRIMARY READS
./runall.sh





##################################################
##################################################
##################################################
# Read segregation
# NOTE: SAM files must be unzipped.

mkdir "${PROJDIR}"/segregation
cd "${PROJDIR}"/segregation

##################################################
# Run segregation
./runall.sh

##################################################
# 160808
# Get ref-both and alt-both counts
# I.e. Get valid pairs w/ at least one end being ambiguously allele-specific
./runall.refalt.sh





##################################################
##################################################
##################################################
# Generate heatmaps

mkdir "${PROJDIR}"/heatmap
cd "${PROJDIR}"/heatmap

#####
# Ensure interaction counts are zipped
find "${DATADIR}"/interactionCounts/ -type f -print0 | xargs -0 gzip -f &


##################################################
# Generate heat maps
./runall.sh	