#!/bin/bash -x
#################################################
#
# This script coordinates 
# 1. the segregation reads
# 2. their pairing
# 3. their binning
# It submits jobs to cluster
# This script copies all files to local /tmp and process them in local storage. 
#
#################################################

set -o nounset
set -o pipefail
set -o errexit 

DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"
BINDIR="<binariesDirName"

RAWDIR=$DATADIR/fullReads
MAPPEDDIR=$DATADIR/samPrimaryAndUnmappedReads
SORTEDMAPPEDDIR=$DATADIR/sortedMappedReads
EXTRACTEDPAIRDIR=$DATADIR/extractedPairs
CLEANPAIRDIR=$DATADIR/cleanedPairs
MIDPOINTS=${DATADIR}/midPoints
FRAGMENTMAPPABILITY=${DATADIR}/fragmentMappability
INTERACTIONCOUNTS=${DATADIR}/interactionCounts

export PATH=$PATH:"${BINDIR}"/UCSC.utilities
assembly=mm10
if [[ ! -e $assembly.size ]]; then
	fetchChromSizes $assembly > $assembly.sizes
fi

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs)) # libIDs to array
IFS=$old_IFS

# mm10 validated snps
snpFileName=spretus.patski.snps.vcf.gz
snpFile=${REFDIR}/mm10pseudoSpretus/${snpFileName}

# Steps to run
runSeg=1
runSAM2mapped=1
runEither=1
runCombine=1
runClean=1
runBinning=1
runSummaries=1

for (( i=0; i<${#libIDs[@]}; i++ )); do
	libID=${libIDs[$i]}
	
	root=$MAPPEDDIR/$libID
	mkdir -p $root
	segroot=$MAPPEDDIR/$libID
	mkdir -p $segroot
	sortedroot=$SORTEDMAPPEDDIR/$libID
	mkdir -p $sortedroot
	extractedroot=$EXTRACTEDPAIRDIR/$libID
	mkdir -p $extractedroot
	cleanedroot=$CLEANPAIRDIR/$libID
	mkdir -p $cleanedroot
	fragmaproot=$FRAGMENTMAPPABILITY/$libID
	mkdir -p $fragmaproot
	intercountroot=$INTERACTIONCOUNTS/$libID
	mkdir -p $intercountroot
	mkdir -p "${MIDPOINTS}"
	
	laneNo=1
	for l in `ls $RAWDIR/$libID/*_1.fq.gz`; do
		laneName="L$laneNo"
		for e in 1 2; do
			base="${laneName}_${e}"
			
			# Segregation step
			if [[ "${runSeg}" -eq 1 ]]; then
			
				jobfile=$libID.$base.segregation.job
				echo hostname > $jobfile
				echo source /etc/profile.d/modules.sh >> $jobfile
				echo module load modules modules-init modules-gs modules-noble >> $jobfile
				echo module load bwa/0.7.3 samtools/1.3 >> $jobfile
				echo module load python/2.7.2 >> $jobfile
				
				# copy black6 sam
				samFileName=$base.sam.gz
				samFile=$MAPPEDDIR/$libID.black6/$samFileName
				inFileName=$base.black6.sam.gz
				echo "if [[ \${TMPDIR:-} != \"\" ]]; then" >> $jobfile
				echo "if [[ ! -s \${TMPDIR:-}/$inFileName ]]; then" >> $jobfile
				echo date >> $jobfile
				echo echo starting copy $samFile >> $jobfile
				echo cp $samFile \${TMPDIR:-}/$inFileName >> $jobfile
				echo date >> $jobfile
				echo echo finished copy >> $jobfile
				echo "fi" >> $jobfile
				echo in1=\${TMPDIR:-}/$inFileName >> $jobfile
				echo "fi" >> $jobfile
				
				# copy spretus sam
				samFileName=$base.sam.gz
				samFile=$MAPPEDDIR/$libID.spretus/$samFileName
				inFileName=$base.spretus.sam.gz
				echo "if [[ \${TMPDIR:-} != \"\" ]]; then" >> $jobfile
				echo "if [[ ! -s \${TMPDIR:-}/$inFileName ]]; then" >> $jobfile
				echo date >> $jobfile
				echo echo starting copy $samFile >> $jobfile
				echo cp $samFile \${TMPDIR:-}/$inFileName >> $jobfile
				echo date >> $jobfile
				echo echo finished copy >> $jobfile
				echo "fi" >> $jobfile
				echo in2=\${TMPDIR:-}/$inFileName >> $jobfile
				echo "fi" >> $jobfile
				
				# copy SNPs
				echo "if [[ \${TMPDIR:-} != \"\" ]]; then" >> $jobfile
				echo "if [[ ! -s \${TMPDIR:-}/$snpFileName ]]; then" >> $jobfile
				echo date >> $jobfile
				echo echo starting copy $snpFile >> $jobfile
				echo cp $snpFile \${TMPDIR:-} >> $jobfile
				echo date >> $jobfile
				echo echo finished copy >> $jobfile
				echo "fi" >> $jobfile
				echo in3=\${TMPDIR:-}/$snpFileName >> $jobfile
				echo "fi" >> $jobfile
				
				# segregation
				echo date >> $jobfile
				echo echo starting segregation>> $jobfile
				echo "./segregate-mapped-reads-py \$in1 \$in2 \$in3 \${TMPDIR:-}/$base > $segroot/$base.txt" >> $jobfile
				echo date >> $jobfile
				echo echo finished segregation >> $jobfile
				
				# copy results files
				echo "if [[ \${TMPDIR:-} != \"\" ]]; then" >> $jobfile
				echo date >> $jobfile
				echo echo starting copy results >> $jobfile
				# for outFileType in alt ref both-ref both-alt; do
				for outFileType in alt ref both-ref; do
				echo out1=\${TMPDIR:-}/$base.$outFileType.sam >> $jobfile
				echo gzip -f \$out1 >> $jobfile
				echo cp \$out1.gz $segroot/ >> $jobfile
				done
				echo date >> $jobfile
				echo echo finished copy results >> $jobfile
				echo fi >> $jobfile
				
				# submit job
				qsub -l mfree=8G -cwd -j y $jobfile
			fi # runSeg
			
			# This step reorders mapped reads
			if [[ "${runSAM2mapped}" -eq 1 ]]; then
				mkdir -p $sortedroot
				for genome in ref alt both-ref; do
					samfile=$segroot/$base.$genome.sam.gz
					outfile=$sortedroot/$base.$genome.gz
					
					jobfile=$libID.$base.$genome.sam2mapped.job
					echo hostname > $jobfile
					echo source /etc/profile.d/modules.sh >> $jobfile
					echo module load modules modules-init modules-gs modules-noble >> $jobfile
					echo module load samtools/0.1.18 >> $jobfile
					echo "samtools view -S $samfile | awk -v genome=$genome -f sam2mapped.awk | sort -k1,1 | gzip > $outfile" >> $jobfile
					
					# submit job
					qsub -l mfree=4G -hold_jid $libID.$base.segregation.job -cwd -j y $jobfile
				done
			 fi # runSAM2mapped

			# This step recombines the segregated reads
			if [[ "${runEither}" -eq 1 ]]; then
				jobfile=$libID.$base.either.job
				echo hostname > $jobfile
				echo "zcat $sortedroot/$base.ref.gz $sortedroot/$base.alt.gz $sortedroot/$base.both-ref.gz | sort -k1,1 | gzip > $sortedroot/$base.either.gz" >> $jobfile
				
				# submit job
				qsub -l mfree=4G -hold_jid $libID.$base.ref.sam2mapped.job,$libID.$base.alt.sam2mapped.job,$libID.$base.both-ref.sam2mapped.job -cwd -j y $jobfile
			fi # runEither
			
		done # for e in 1 2; do
		
		# This step pairs up read ends
		if [[ "${runCombine}" -eq 1 ]]; then

			file1=$sortedroot/${laneName}_1.either.gz
			file2=$sortedroot/${laneName}_2.either.gz
			outfile=$extractedroot/${laneName}
		
			jobfile=${libID}.$laneName.combine.job
			echo "#!/bin/bash -ex" > $jobfile
			echo "#$ -cwd" >> $jobfile
			echo "#$ -j y" >> $jobfile
			echo hostname >> $jobfile
			echo "join <(zcat $file1)	<(zcat $file2) > $outfile" >> $jobfile
			echo "n=\`cat $outfile | wc -l \`" >> $jobfile
			echo "echo -e \"${libID}\t${laneName}\t\$n \" " >> $jobfile
			echo "gzip -f $outfile" >> $jobfile
			qsub -l mfree=4G -hold_jid $libID.${laneName}_1.either.job,$libID.${laneName}_2.either.job -cwd -j y $jobfile
		fi # runCombine
		
		# This step cleans up paired reads, sorting the columns
		if [[ "${runClean}" -eq 1 ]]; then
			infile=$extractedroot/${laneName}.gz
			outfile=$cleanedroot/${laneName}
		
			jobfile=${libID}.$laneName.clean.dup.job
			afterjobfile=${libID}.clean.dup.out
			echo "#!/bin/bash -ex" > $jobfile
			echo "#$ -cwd" >> $jobfile
			echo "#$ -l mfree=8G" >> $jobfile
			echo "#$ -j y" >> $jobfile
			echo source /etc/profile.d/modules.sh >> $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo module load python/2.7.2 >> $jobfile
			echo hostname >> $jobfile
			echo "touch $outfile" >> $jobfile
			echo "" >> $jobfile
			echo "zcat $infile | python sort_loci.py - 2 | awk '{s1=\"+\"; s2=\"+\"; if (\$2==\"16\") {s1=\"-\"}; if (\$7==\"16\") {s2=\"-\"}; printf(\"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\", \$3, \$4, s1, \$8, \$9, s2, \$5, \$6, \$10, \$11, \$1)}' | sort -u -k1,6 | awk 'BEGIN{OFS=\"\\t\";}{ print \$11,\$3,\$1,\$2,\$7,\$8,\$6,\$4,\$5,\$9,\$10;}' >> $outfile" >> $jobfile
			echo "gzip -f $outfile " >> $jobfile
			
			# submit job
			qsub -l mfree=4G -hold_jid $libID.$laneName.combine.job -cwd -j y $jobfile
		fi # runClean
		
		laneNo=$(($laneNo+1))
		
	done # for l in `ls $RAWDIR/$libID/*_1.fq.gz`;
	
	# Now bin the reads
	# all lanes will be joined, so no need for looping for all lanes
	if [[ "${runBinning}" -eq 1 ]]; then
		for BINSIZE in 500000 100000 40000; do
			jobfile=${libID}.${BINSIZE}.interation.count.job
			echo "#!/bin/bash -ex" > $jobfile
			echo "#$ -cwd" >> $jobfile
			echo "#$ -j y" >> $jobfile
			echo source /etc/profile.d/modules.sh >> $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo module load python/2.7.2 >> $jobfile
			echo hostname >> $jobfile
			
			output_file=${MIDPOINTS}/$assembly.midPoints.${BINSIZE}
			outfile1=$intercountroot/${BINSIZE}
			outfile2=$fragmaproot/${BINSIZE}
			echo "python generate_binned_midpoints.py $BINSIZE $assembly.sizes $output_file" >> $jobfile
			
			infile1=$cleanedroot/*.gz
			infile2=$output_file
			MAPPABILITYTHRESHOLD=1
			DISTANCETHRESHOLD=20000
			echo "python count_interactions_per_binned_fragPairs.py <(zcat ${infile1}) $infile2 $BINSIZE $MAPPABILITYTHRESHOLD $DISTANCETHRESHOLD $outfile1 $outfile2 > $libID.$BINSIZE.summary.txt" >> $jobfile
			echo "cat ${outfile1}.either | sort -nr -k7,7 | sort -s -k1,4 > $outfile1.either">> $jobfile
			
			# submit job
			JobIDs="$libID.*.clean.dup.job"
			echo $JobIDs
			qsub -N $jobfile -l mfree=4G -hold_jid "$JobIDs" -cwd -j y $jobfile
			
		done # for BINSIZE
	fi # runBinning
	
done # for (( i=0; i<${#libIDs[@]}; i++ )); do


if [[ "${runSummaries}" -eq 1 ]]; then
	jobfile=interaction.stats.job
	echo hostname > $jobfile
	echo export PATH=\$PATH:$BINDIR >> $jobfile
	echo "./summarize-interaction-stats.sh cleanedPairs > interaction.stats.txt" >> $jobfile
	echo "rdb2html -noformatline interaction.stats.txt > interaction.stats.html" >> $jobfile
	echo "./summarize-interaction-segregation-stats.sh cleanedPairs > interaction.segregation.stats.txt" >> $jobfile
	echo "cat interaction.segregation.stats.txt | sed 's/both-ref/both/g' > interaction.segregation.stats.tmp"	>> $jobfile
	echo "mv interaction.segregation.stats.tmp interaction.segregation.stats.txt"	>> $jobfile
	echo "rdb2html -noformatline interaction.segregation.stats.txt > interaction.segregation.stats.html" >> $jobfile
	
	# submit job
	JobIDs="$libID.*.clean.dup.job"
	echo $JobIDs
	qsub -l mfree=32G -l h_rt=47:59:59 -hold_jid "$JobIDs" -cwd -j y $jobfile
fi # runSummaries
