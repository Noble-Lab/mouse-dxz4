#!/bin/bash -ex
#################################################
#
# This script generates a job script to map reads using bwa mem and submits it to cluster.
# This script copies all files to local /tmp and process them in local storage. 
#
#################################################

LIB_ID=$1 
INDICES=$2
ASSEMBLY=$3
DATADIR=$4

RAWDIR=$DATADIR/fullReads/$1
MAPPEDDIR=$DATADIR/mappedReads/$1.$ASSEMBLY

mkdir -p $MAPPEDDIR

laneNo=1
for i in $(ls $RAWDIR/*_1.fq.gz); do
	echo $i
	filename=${i##*/}
	basename=${filename%%_1.fq.gz}
	echo $basename
	laneName="L$laneNo"
	for e in 1 2; do
		base="${laneName}_${e}"
		
		jobfile=${LIB_ID}.$base.$ASSEMBLY.mm10.step1.job
		
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
		echo module load bwa/0.7.3 samtools/1.3 >> $jobfile
		echo "LC_COLLATE=C; LC_ALL=C ; LANG=C ; export LC_ALL LANG">> $jobfile
		echo "echo \"Environment: LC_COLLATE=\$LC_COLLATE, LC_ALL = \$LC_ALL, LANG = \$LANG \" " >> $jobfile
		echo "" >> $jobfile
		echo hostname >> $jobfile
		echo "" >> $jobfile
			
		# If on cluster, copy data files to temporary directory.
		echo "if [[ \${TMPDIR:-} != \"\" ]]; then" >> $jobfile
			echo date >> $jobfile
			echo echo starting copy >> $jobfile
			echo cp $fastqFile \${TMPDIR:-} >> $jobfile
			echo date >> $jobfile
			echo echo finished copy >> $jobfile
			echo in1=\${TMPDIR:-}/$fastqFileName >> $jobfile
			echo out1=\${TMPDIR:-}/$samFileName >> $jobfile
		echo else >> $jobfile
			echo in1=$fastqFile >> $jobfile
			echo out1=$samFile >> $jobfile
		echo fi >> $jobfile
		
		# mapping
		echo date >> $jobfile
		echo echo starting mapping $fastqFile >> $jobfile
		echo "bwa mem -t 8 -M $INDICES <(zcat \$in1) > \$out1" >> $jobfile
		echo date >> $jobfile
		echo echo finished mapping >> $jobfile
		# copy results files
		echo "if [ \${TMPDIR:-} != \"\" ]; then" >> $jobfile
			echo date >> $jobfile
			echo echo starting copy results >> $jobfile
			echo cp \$out1 $MAPPEDDIR/ >> $jobfile
			echo date >> $jobfile
			echo echo finished copy results >> $jobfile
		echo fi >> $jobfile
		
		# Submit job
		qsub -pe serial 8 -l mfree=4G -q noble-long.q -cwd -j y -N $jobfile $jobfile
	
	done
	
	laneNo=$(($laneNo+1))
done

