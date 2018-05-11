#!/bin/bash -ex
#################################################
#
# This is script copy all files to local /tmp and process them in local storage. 
# It won't copy *.sam file. Instead, will copy *.sai and generate it locally.
#
# 160414 Giancarlo: Made code slightly cleaner
#
# 160418 Giancarlo: Parse bwa mem results
#
#################################################

LIB_ID=$1 
ASSEMBLY=$2
DATADIR=$3

RAWDIR=$DATADIR/fullReads/$1
MAPPEDDIR=$DATADIR/mappedReadsPE/$1.$ASSEMBLY

SORTEDDIR=$DATADIR/samPrimaryAndUnmappedReadsPE/$1.$ASSEMBLY
mkdir -p $SORTEDDIR

laneNo=1
for i in `ls $RAWDIR/*_1.fq.gz`; do
	echo $i
	filename=${i##*/}
	basename=${filename%%_1.fq.gz}
	echo $basename
	laneName="L$laneNo"
	base="${laneName}"
		jobfile=${LIB_ID}.$base.$ASSEMBLY.step2c.job
	
	samFileName=$base.bam
	samFile=$MAPPEDDIR/$samFileName

	echo "#!/bin/bash" > $jobfile
	echo "#$ -cwd" >> $jobfile
	echo "#$ -l mfree=16G" >> $jobfile
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
	echo "samtools view -F 256 \$in1 | sort -k 1,1 | gzip -f > \$out1" >> $jobfile
	echo date >> $jobfile
	echo "echo finished filtering mapped reads" >> $jobfile
	
	echo >> $jobfile
	
	###########################################################
	## Submit job
	qsub -l h_rt=7:59:59 -hold_jid ${LIB_ID}.$base.$ASSEMBLY.step1.job -N $jobfile $jobfile -cwd -j y 

laneNo=$(($laneNo+1))
done

