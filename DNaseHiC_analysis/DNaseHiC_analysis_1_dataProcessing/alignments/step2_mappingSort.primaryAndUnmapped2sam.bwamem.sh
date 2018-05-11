#!/bin/bash -ex
#################################################
#
# This script filters out secondary aligned reads and sorts them.
#
#################################################

LIB_ID=$1 
ASSEMBLY=$2
DATADIR=$3

RAWDIR=$DATADIR/fullReads/$1
MAPPEDDIR=$DATADIR/mappedReads/$1.$ASSEMBLY

SORTEDDIR=$DATADIR/samPrimaryAndUnmappedReads/$1.$ASSEMBLY
mkdir -p $SORTEDDIR

laneNo=1
for i in `ls $RAWDIR/*_1.fq.gz`; do
	echo $i
	filename=${i##*/}
	basename=${filename%%_1.fq.gz}
	echo $basename
	laneName="L$laneNo"
	for e in 1 2; do
		base="${laneName}_${e}"

		jobfile=${LIB_ID}.$base.$ASSEMBLY.step2c.job
		
		samFileName=$base.bam
		samFile=$MAPPEDDIR/$samFileName
		fastqFileName=${basename}_${e}.fq.gz
		fastqFile=$RAWDIR/$fastqFileName

		echo "#!/bin/bash" > $jobfile
		echo "#$ -cwd" >> $jobfile
		echo "#$ -l h_vmem=16G" >> $jobfile
		echo "#$ -j y" >> $jobfile
		echo source /etc/profile.d/modules.sh >> $jobfile
		echo module load modules modules-init modules-gs modules-noble >> $jobfile
		echo module load samtools/1.3 >> $jobfile
		echo "LC_COLLATE=C; LC_ALL=C ; LANG=C ; export LC_ALL LANG">> $jobfile
		echo "echo \"Environment: LC_COLLATE=\$LC_COLLATE, LC_ALL = \$LC_ALL, LANG = \$LANG \" " >> $jobfile
		echo "" >> $jobfile
		echo hostname >> $jobfile
		echo "" >> $jobfile
	
		# Drop seconday reads and sort
		outFileName=$base.sam.gz
		outFile=$SORTEDDIR/$outFileName
		echo in1=$samFile >> $jobfile	
		echo out1=$outFile >> $jobfile
		echo date >> $jobfile
		echo "echo starting filtering for mapped and unmapped" >> $jobfile
		echo "samtools view -h -F 256 \$in1 | sort -k 1,1 | gzip -f > \$out1" >> $jobfile
		echo date >> $jobfile
		echo "echo finished filtering mapped reads" >> $jobfile
		echo >> $jobfile

		# Submit job
		qsub -l h_rt=7:59:59 -hold_jid ${LIB_ID}.$base.$ASSEMBLY.step1.job -N $jobfile $jobfile -cwd -j y 
	
	done
	
	laneNo=$(($laneNo+1))
done

