#!/bin/bash -x
#################################################
#
# 160809 Giancarlo
# Added calls to parse, ICE, and plot the ref-ref/both and alt-alt/both read pair combinations matrices
#
# 170111 Giancarlo
# mm10
#
#################################################

set -o nounset
set -o pipefail
set -o errexit 

DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"
BINDIR="<binariesDirName"

INTERACTIONCOUNTSDIR=$DATADIR/interactionCounts
CLEANEDPAIRSDIR=$DATADIR/cleanedPairs

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble

module load lmutil/10.8.5
module load matlab/2008b

module load python/2.7.2 
module load samtools/0.1.18
module load bowtie2/2.0.0-beta6
module load numpy/1.6.1 
module load scipy/0.10.0
module load zlib/1.2.6
module load hdf5/1.8.3

export PATH=$PATH:$BINDIR
export PYTHONPATH=$PYTHONPATH:$BINDIR/mirnylab-mirnylib-213b592a302c:$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/src

# GB 170111
assembly=mm10
if [ ! -e mm10.sizes ]; then
	fetchChromSizes mm10 > mm10.sizes
fi

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs)) # libIDs to array
IFS=$old_IFS

resolutions=(500000 10000 40000)

RUNSTEP1=1
RUNSTEP2=1

# Parse Hi-C data and generate binned matrices in hdf5 format
if [[ "${RUNSTEP1}" -eq 1 ]]; then
	for libID in ${libIDs[@]}; do
		genomePath=$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/fasta/$assembly
		for resolution in ${resolutions[@]}; do
	
			dir=$libID.$resolution
			mkdir -p $dir
		
			if [[ "${resolution}" -lt 100000 ]]; then
					VMEM="64G"
			elif [[ "${resolution}" -lt 500000 ]]; then
					VMEM="32G"
			else 
					VMEM="16G"
			fi
			
			jobfile=$libID.$assembly.$resolution.parse.job
			echo source /etc/profile.d/modules.sh > $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo module load lmutil/10.8.5 >> $jobfile
			echo module load matlab/2008b >> $jobfile
			echo module load python/2.7.2 >> $jobfile
			echo module load samtools/0.1.18 >> $jobfile
			echo module load bowtie2/2.0.0-beta6 >> $jobfile
			echo module load numpy/1.6.1 >> $jobfile
			echo module load scipy/0.10.0 >> $jobfile
			echo module load zlib/1.2.6 >> $jobfile
			echo module load hdf5/1.8.3 >> $jobfile
			# GB 170111
			echo export PYTHONPATH=\$PYTHONPATH:\$BINDIR/mirnylab-mirnylib-213b592a302c:\$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/src >> $jobfile
			
			##############################################
			# Either		
			echo ./parsehic.either.py $assembly.sizes M $INTERACTIONCOUNTSDIR/$libID/$resolution.either.gz $dir $genomePath $resolution >> $jobfile
			
			##############################################
			# 160809 
			# ref-both/ref and alt-both/alt count matrices
			echo ./parsehic.either.py $assembly.sizes M $INTERACTIONCOUNTSDIR/$libID/$resolution.ref.gz $dir $genomePath $resolution >> $jobfile
			echo ./parsehic.either.py $assembly.sizes M $INTERACTIONCOUNTSDIR/$libID/$resolution.alt.gz $dir $genomePath $resolution >> $jobfile

			##############################################
			# All other end combinations
			for genome1 in ref alt both-ref; do
				for genome2 in ref alt both-ref; do
					if [[ ! $genome1 > $genome2 ]]; then
						echo ./parsehic.genome.py $assembly.sizes M $INTERACTIONCOUNTSDIR/$libID/$resolution.gz $dir $genomePath $resolution $genome1 $genome2 >> $jobfile
						#fi
					fi
				done
			done
		
			# submit job
			qsub -l mfree=${VMEM} -cwd -j y $jobfile
		
		done
	done
fi


# Normalization and plots
if [[ "${RUNSTEP2}" -eq 1 ]]; then
	for libID in ${libIDs[@]}; do
		for resolution in ${resolutions[@]}; do
			dir=$libID.$resolution
		
			if [[ "${resolution}" -lt 100000 ]]; then
					VMEM="64G"
			elif [[ "${resolution}" -lt 500000 ]]; then
					VMEM="32G"
			else 
					VMEM="16G"
			fi
		
			##############################################
			# Either		
			jobfile=$libID.$assembly.$resolution.either.IC.job
			echo source /etc/profile.d/modules.sh > $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo module load lmutil/10.8.5 >> $jobfile
			echo module load matlab/2008b >> $jobfile
			echo module load python/2.7.2 >> $jobfile
			echo module load samtools/0.1.18 >> $jobfile
			echo module load R/2.14.0 >> $jobfile
			echo module load bowtie2/2.0.0-beta6 >> $jobfile
			echo module load numpy/1.6.1 >> $jobfile
			echo module load scipy/0.10.0 >> $jobfile
			echo module load zlib/1.2.6 >> $jobfile
			echo module load hdf5/1.8.3 >> $jobfile
			# GB 170111
			echo export PYTHONPATH=\$PYTHONPATH:\$BINDIR/mirnylab-mirnylib-213b592a302c:\$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/src >> $jobfile
		
			echo ./mirnylab.ICE.heatmap.py $dir/$resolution.either.chr.hdf5 $dir/$resolution.either $assembly $resolution $libID >> $jobfile
			echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.either.chr.hdf5 $dir/$resolution.either $assembly $resolution $libID.raw >> $jobfile
			echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.either.IC.chr.hdf5 $dir/$resolution.either.IC $assembly $resolution $libID.corrected >> $jobfile
		
			# submit job
			qsub -hold_jid $libID.$assembly.$resolution.parse.job -l mfree="${VMEM}" -cwd -j y $jobfile



			##############################################
			# 160809 
			# ref-both/ref and alt-both/alt count matrices
		
			#####
			# ref
			jobfile=$libID.$assembly.$resolution.ref.IC.job
			echo source /etc/profile.d/modules.sh > $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo module load lmutil/10.8.5 >> $jobfile
			echo module load matlab/2008b >> $jobfile
			echo module load python/2.7.2 >> $jobfile
			echo module load samtools/0.1.18 >> $jobfile
			echo module load R/2.14.0 >> $jobfile
			echo module load bowtie2/2.0.0-beta6 >> $jobfile
			echo module load numpy/1.6.1 >> $jobfile
			echo module load scipy/0.10.0 >> $jobfile
			echo module load zlib/1.2.6 >> $jobfile
			echo module load hdf5/1.8.3 >> $jobfile
			# GB 170111
			echo export PYTHONPATH=\$PYTHONPATH:\$BINDIR/mirnylab-mirnylib-213b592a302c:\$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/src >> $jobfile
		
			echo ./mirnylab.ICE.heatmap.py $dir/$resolution.ref.chr.hdf5 $dir/$resolution.ref $assembly $resolution $libID >> $jobfile
			echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.ref.chr.hdf5 $dir/$resolution.ref $assembly $resolution $libID.raw >> $jobfile
			echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.ref.IC.chr.hdf5 $dir/$resolution.ref.IC $assembly $resolution $libID.corrected >> $jobfile
		
			# submit job
			qsub -hold_jid $libID.$assembly.$resolution.parse.job -l mfree="${VMEM}" -cwd -j y $jobfile
		
		 
			#####
			# alt
			jobfile=$libID.$assembly.$resolution.alt.IC.job
			echo source /etc/profile.d/modules.sh > $jobfile
			echo module load modules modules-init modules-gs modules-noble >> $jobfile
			echo module load lmutil/10.8.5 >> $jobfile
			echo module load matlab/2008b >> $jobfile
			echo module load python/2.7.2 >> $jobfile
			echo module load samtools/0.1.18 >> $jobfile
			echo module load R/2.14.0 >> $jobfile
			echo module load bowtie2/2.0.0-beta6 >> $jobfile
			echo module load numpy/1.6.1 >> $jobfile
			echo module load scipy/0.10.0 >> $jobfile
			echo module load zlib/1.2.6 >> $jobfile
			echo module load hdf5/1.8.3 >> $jobfile
			# GB 170111
			echo export PYTHONPATH=\$PYTHONPATH:\$BINDIR/mirnylab-mirnylib-213b592a302c:\$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/src >> $jobfile
		
			echo ./mirnylab.ICE.heatmap.py $dir/$resolution.alt.chr.hdf5 $dir/$resolution.alt $assembly $resolution $libID >> $jobfile
			echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.alt.chr.hdf5 $dir/$resolution.alt $assembly $resolution $libID.raw >> $jobfile
			echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.alt.IC.chr.hdf5 $dir/$resolution.alt.IC $assembly $resolution $libID.corrected >> $jobfile
		
			# submit job
			qsub -hold_jid $libID.$assembly.$resolution.parse.job -l mfree="${VMEM}" -cwd -j y $jobfile
		
		
			##############################################
			# All other end combinations	
			for genome1 in ref alt both-ref; do
				genome1Name=$genome1
				if [ $genome1 == "both-ref" ]; then
					genome1Name="both"
				fi
				for genome2 in ref alt both-ref; do
					genome2Name=$genome2
					if [ $genome2 == "both-ref" ]; then
						genome2Name="both"
					fi
					if [[ ! $genome1 > $genome2 ]]; then
						jobfile=$libID.$assembly.$resolution.$genome1.$genome2.IC.job
						echo source /etc/profile.d/modules.sh > $jobfile
						echo module load modules modules-init modules-gs modules-noble >> $jobfile
						echo module load lmutil/10.8.5 >> $jobfile
						echo module load matlab/2008b >> $jobfile
						echo module load python/2.7.2 >> $jobfile
						echo module load samtools/0.1.18 >> $jobfile
						echo module load R/2.14.0 >> $jobfile
						echo module load bowtie2/2.0.0-beta6 >> $jobfile
						echo module load numpy/1.6.1 >> $jobfile
						echo module load scipy/0.10.0 >> $jobfile
						echo module load zlib/1.2.6 >> $jobfile
						echo module load hdf5/1.8.3 >> $jobfile
						# GB 170111
						echo export PYTHONPATH=\$PYTHONPATH:\$BINDIR/mirnylab-mirnylib-213b592a302c:\$BINDIR/mirnylab-hiclib-9ed8d9e0ca7f/src >> $jobfile
					
						echo ./mirnylab.ICE.byBiases.heatmap.py $dir/$resolution.$genome1.$genome2.chr.hdf5 $dir/$resolution.either.IC.chr.biases.txt $dir/$resolution.$genome1.$genome2 $assembly $resolution $libID >> $jobfile
						echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.$genome1.$genome2.chr.hdf5 $dir/$resolution.$genome1.$genome2 $assembly $resolution $libID.$genome1Name.$genome2Name.raw >> $jobfile
						echo ./mirnylab.HiResHiC.chr.heatmap.py $dir/$resolution.$genome1.$genome2.IC.byBiases.chr.hdf5 $dir/$resolution.$genome1.$genome2.IC.byBiases $assembly $resolution $libID.$genome1Name.$genome2Name.corrected >> $jobfile
					
						# submit job
						qsub -hold_jid $libID.$assembly.$resolution.either.IC.job -l mfree="${VMEM}" -cwd -j y $jobfile
					fi
				done
			done
		done
	done
fi