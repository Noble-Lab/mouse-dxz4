#!/bin/bash -x
######################################################################
#
# This script goes through the existing binned interaction counts file
# and re-tallies ref-both and alt-both counts to get read pairs where
# at least one end is unambiguously aligned to an allele. 
# 
######################################################################

set -o nounset
set -o pipefail
set -o errexit 

export PATH=$PATH:$HOME/bin/UCSC.utilities

DATADIR=$(cat ../DATADIR)

INTERACTIONCOUNTS=${DATADIR}/interactionCounts

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs )) # libIDs to array
IFS=$old_IFS

for (( i=0; i<${#libIDs[@]}; i++ )); do
	libID=${libIDs[$i]}
	
	intercountroot=$INTERACTIONCOUNTS/$libID
	
	for BINSIZE in 500000 100000 40000; do
	
		infile=$intercountroot/${BINSIZE}
		reffile=$intercountroot/${BINSIZE}.ref 
		altfile=$intercountroot/${BINSIZE}.alt
	
		if [[ -s $infile.gz ]]; then
			
			jobfile=${libID}.${BINSIZE}.interation.count.ref.CIGARsed.job
			echo "#!/bin/bash -ex" > $jobfile
			echo "#$ -cwd" >> $jobfile
			echo "#$ -j y" >> $jobfile
			echo source /etc/profile.d/modules.sh >> $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo PATH=\$PATH:/net/noble/vol1/home/noble/bin >> $jobfile

			echo hostname >> $jobfile
			
			echo "zcat $infile.gz | \\" >> $jobfile
			echo "awk '{ \\" >> $jobfile
			echo "locusID=sprintf(\"%s-%s-%s-%s\", \$1, \$2, \$3, \$4); \\" >> $jobfile
			echo "if((\$5~/ref/ || \$5~/both/) && (\$6~/ref/ || \$6~/both/) && !(\$5~/both/ && \$6~/both/)){refhash[locusID]+=\$7} \\" >> $jobfile
			echo "} \\" >> $jobfile
			echo "END{ \\" >> $jobfile
			echo "for(key in refhash){split(key, keyArray, \"-\"); printf(\"%s\t%s\t%s\t%s\t%s\n\", keyArray[1], keyArray[2], keyArray[3], keyArray[4], refhash[key]) } \\" >> $jobfile
			echo "}' | gzip -f > ${reffile}.gz" >> $jobfile
			
			chmod u+g $jobfile
			qsub -l h_vmem=16G -l h_rt="7:59:59" -cwd $jobfile
			
						
			jobfile=${libID}.${BINSIZE}.interation.count.alt.CIGARsed.job
			echo "#!/bin/bash -ex" > $jobfile
			echo "#$ -cwd" >> $jobfile
			echo "#$ -j y" >> $jobfile
			echo source /etc/profile.d/modules.sh >> $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo hostname >> $jobfile
			
			echo "zcat $infile.gz | \\" >> $jobfile
			echo "awk '{ \\" >> $jobfile
			echo "locusID=sprintf(\"%s-%s-%s-%s\", \$1, \$2, \$3, \$4); \\" >> $jobfile
			echo "if((\$5~/alt/ || \$5~/both/) && (\$6~/alt/ || \$6~/both/) && !(\$5~/both/ && \$6~/both/)){althash[locusID]+=\$7} \\" >> $jobfile
			echo "} \\" >> $jobfile
			echo "END{ \\" >> $jobfile
			echo "for(key in althash){split(key, keyArray, \"-\"); printf(\"%s\t%s\t%s\t%s\t%s\n\", keyArray[1], keyArray[2], keyArray[3], keyArray[4], althash[key]) } \\" >> $jobfile
			echo "}' | gzip -f > ${altfile}.gz" >> $jobfile
			
			# submit job
			chmod u+g $jobfile
			qsub -l h_vmem=4G -l h_rt="7:59:59" -cwd $jobfile

		else 
			echo "$infile.gz does not exist."
		fi 

	done
done
